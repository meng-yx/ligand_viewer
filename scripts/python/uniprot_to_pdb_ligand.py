#!/usr/bin/env python3
"""Build UniProt -> PDB ligand mapping and stream results to CSV."""

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import requests
from joblib import Parallel, delayed
from tqdm import tqdm

RCSB_CORE_API = "https://data.rcsb.org/rest/v1/core"

OUTPUT_COLUMNS = ["uniprot_id", "pdb_id", "ligand_code", "ligand_name", "smiles"]


def read_uniprot_ids(path: Path) -> List[str]:
    """Read UniProt IDs from a text file and deduplicate by first appearance."""
    ordered_ids: List[str] = []
    seen: Set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            value = value.upper()
            if value in seen:
                continue
            seen.add(value)
            ordered_ids.append(value)
    return ordered_ids


def read_ligand_blacklist(path: Optional[Path]) -> Set[str]:
    """Read ligand codes to ignore from text file.

    Supports either:
      - "CODE"
      - "CODE, ligand name"
    """
    if path is None:
        return set()
    codes: Set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            if value.startswith("#"):
                continue
            # Keep only ligand code (first CSV field) so lines like
            # "MG, MAGNESIUM ION" blacklist code "MG".
            code = value.split(",", 1)[0].strip().upper()
            if code:
                codes.add(code)
    return codes


def get_pdb_accessions(uniprot_id: str) -> List[str]:
    """Return unique PDB IDs cross-referenced for a UniProt accession."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
    except (requests.RequestException, ValueError) as exc:
        print(f"{uniprot_id}: UniProt fetch failed ({exc})")
        return []

    pdb_ids: List[str] = []
    seen: Set[str] = set()
    for dbref in data.get("uniProtKBCrossReferences", []):
        if dbref.get("database") != "PDB":
            continue
        pdb_id = dbref.get("id")
        if pdb_id and pdb_id not in seen:
            seen.add(pdb_id)
            pdb_ids.append(pdb_id.upper())
    return pdb_ids


def _get_chemcomp_smiles(comp_id: str) -> Optional[str]:
    """Fetch SMILES for a chemical component ID."""
    try:
        response = requests.get(f"{RCSB_CORE_API}/chemcomp/{comp_id}", timeout=30)
        response.raise_for_status()
        descriptor = response.json().get("rcsb_chem_comp_descriptor") or {}
        return descriptor.get("SMILES_stereo") or descriptor.get("SMILES")
    except requests.RequestException:
        return None


def _build_polymer_uniprot_maps(
    pdb_id: str, polymer_entity_ids: Sequence[str]
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Map polymer entity_id and asym_id -> UniProt accessions for a PDB entry."""
    entity_to_uniprot: Dict[str, List[str]] = {}
    asym_to_uniprot: Dict[str, List[str]] = {}

    for entity_id in polymer_entity_ids:
        try:
            response = requests.get(
                f"{RCSB_CORE_API}/polymer_entity/{pdb_id}/{entity_id}",
                timeout=30,
            )
            if response.status_code != 200:
                continue
            container = response.json().get("rcsb_polymer_entity_container_identifiers") or {}
        except (requests.RequestException, ValueError, TypeError, KeyError):
            continue
        uniprot_ids = container.get("uniprot_ids") or []
        entity_to_uniprot[str(entity_id)] = uniprot_ids

        for asym_id in (container.get("asym_ids") or []) + (
            container.get("auth_asym_ids") or []
        ):
            asym_to_uniprot[asym_id] = uniprot_ids

    return entity_to_uniprot, asym_to_uniprot


def _neighbor_uniprot_ids(
    neighbor: Dict[str, str],
    entity_to_uniprot: Dict[str, List[str]],
    asym_to_uniprot: Dict[str, List[str]],
) -> Set[str]:
    """Resolve a target neighbor record to UniProt accessions."""
    uniprot_ids: Set[str] = set()

    target_entity_id = neighbor.get("target_entity_id")
    if target_entity_id is not None:
        for uid in entity_to_uniprot.get(str(target_entity_id), []):
            uniprot_ids.add(uid)

    target_asym_id = neighbor.get("target_asym_id")
    if target_asym_id is not None:
        for uid in asym_to_uniprot.get(target_asym_id, []):
            uniprot_ids.add(uid)

    return uniprot_ids


def _is_subject_of_investigation(entity_data: Dict) -> bool:
    """True if the non-polymer entity is annotated as SUBJECT_OF_INVESTIGATION."""
    for annotation in entity_data.get("rcsb_nonpolymer_entity_annotation") or []:
        if annotation.get("type") == "SUBJECT_OF_INVESTIGATION":
            return True

    for feature in entity_data.get("rcsb_nonpolymer_entity_feature") or []:
        if feature.get("type") == "SUBJECT_OF_INVESTIGATION":
            return True

    return False


def get_ligands_from_pdb(
    pdb_id: str, uniprot_id: str, filter_by_subject_of_investigation: bool = False
) -> List[Dict[str, Optional[str]]]:
    """
    Query RCSB PDB API for ligands in a structure that contact the query protein.

    Returns a list of dicts with: ligand_code, ligand_name, smiles.
    """
    pdb_id = pdb_id.upper()
    uniprot_id = uniprot_id.upper()

    try:
        entry_response = requests.get(f"{RCSB_CORE_API}/entry/{pdb_id}", timeout=30)
        if entry_response.status_code != 200:
            return []
        container = entry_response.json().get("rcsb_entry_container_identifiers") or {}
    except (requests.RequestException, ValueError, TypeError, KeyError):
        return []

    nonpolymer_entity_ids = container.get("non_polymer_entity_ids") or []
    polymer_entity_ids = container.get("polymer_entity_ids") or []
    if not nonpolymer_entity_ids:
        return []

    entity_to_uniprot, asym_to_uniprot = _build_polymer_uniprot_maps(
        pdb_id, polymer_entity_ids
    )

    ligand_rows: List[Dict[str, Optional[str]]] = []
    seen_comp_ids: Set[str] = set()

    for entity_id in nonpolymer_entity_ids:
        try:
            entity_response = requests.get(
                f"{RCSB_CORE_API}/nonpolymer_entity/{pdb_id}/{entity_id}",
                timeout=30,
            )
            if entity_response.status_code != 200:
                continue
            entity_data = entity_response.json()
        except (requests.RequestException, ValueError, TypeError, KeyError):
            continue
        comp_id = (entity_data.get("pdbx_entity_nonpoly") or {}).get("comp_id")
        comp_name = (entity_data.get("pdbx_entity_nonpoly") or {}).get("name")
        if not comp_id or comp_id in seen_comp_ids:
            continue

        if filter_by_subject_of_investigation and not _is_subject_of_investigation(entity_data):
            continue

        instance_asym_ids = (
            (entity_data.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get(
                "asym_ids"
            )
            or []
        )

        binds_query_protein = False
        for asym_id in instance_asym_ids:
            try:
                instance_response = requests.get(
                    f"{RCSB_CORE_API}/nonpolymer_entity_instance/{pdb_id}/{asym_id}",
                    timeout=30,
                )
                if instance_response.status_code != 200:
                    continue
                neighbors = instance_response.json().get("rcsb_target_neighbors") or []
            except (requests.RequestException, ValueError, TypeError, KeyError):
                continue
            for neighbor in neighbors:
                neighbor_uniprots = _neighbor_uniprot_ids(
                    neighbor, entity_to_uniprot, asym_to_uniprot
                )
                if uniprot_id in neighbor_uniprots:
                    binds_query_protein = True
                    break
            if binds_query_protein:
                break

        if not binds_query_protein:
            continue

        seen_comp_ids.add(comp_id)
        ligand_rows.append(
            {
                "ligand_code": comp_id,
                "ligand_name": comp_name,
                "smiles": _get_chemcomp_smiles(comp_id),
            }
        )

    return ligand_rows


def build_rows_for_uniprot(
    uniprot_id: str, ligand_blacklist: Set[str]
) -> List[Dict[str, Optional[str]]]:
    """Build output rows for one UniProt accession."""
    output_rows: List[Dict[str, Optional[str]]] = []
    pdb_ids = get_pdb_accessions(uniprot_id)
    for pdb_id in pdb_ids:
        try:
            ligands = get_ligands_from_pdb(pdb_id, uniprot_id)
        except (requests.RequestException, TypeError, KeyError, AttributeError) as exc:
            print(f"{uniprot_id} {pdb_id}: skipped ({exc})")
            continue

        filtered_ligands = [
            ligand
            for ligand in ligands
            if (ligand.get("ligand_code") or "").upper() not in ligand_blacklist
        ]

        if not filtered_ligands:
            output_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "pdb_id": pdb_id,
                    "ligand_code": None,
                    "ligand_name": None,
                    "smiles": None,
                }
            )
            continue

        for ligand in filtered_ligands:
            output_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "pdb_id": pdb_id,
                    "ligand_code": ligand.get("ligand_code"),
                    "ligand_name": ligand.get("ligand_name"),
                    "smiles": ligand.get("smiles"),
                }
            )
    return output_rows


def iter_parallel_results(
    uniprot_ids: Sequence[str], ligand_blacklist: Set[str], n_jobs: int
) -> Iterable[List[Dict[str, Optional[str]]]]:
    """Yield per-UniProt result rows as tasks complete."""
    parallel = Parallel(n_jobs=n_jobs, prefer="threads", return_as="generator_unordered")
    generator = parallel(
        delayed(build_rows_for_uniprot)(uniprot_id, ligand_blacklist)
        for uniprot_id in uniprot_ids
    )
    yield from generator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build UniProt -> PDB ligand map and stream rows to CSV."
    )
    parser.add_argument("--input", required=True, help="Input .txt with one UniProt ID per row.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument(
        "--blacklist",
        default=None,
        help="Optional .txt containing ligand codes to exclude (one per row).",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="joblib worker count (default: -1 uses all cores).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    blacklist_path = Path(args.blacklist) if args.blacklist else None

    print(f"Input path: {input_path}")
    print(f"Output path: {output_path}")
    print(f"Blacklist path: {blacklist_path}")

    uniprot_ids = read_uniprot_ids(input_path)
    ligand_blacklist = read_ligand_blacklist(blacklist_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    written_rows = 0

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for rows in tqdm(
            iter_parallel_results(uniprot_ids, ligand_blacklist, args.n_jobs),
            total=len(uniprot_ids),
            desc="Proteins",
            unit="protein",
        ):
            if not rows:
                continue
            writer.writerows(rows)
            handle.flush()
            written_rows += len(rows)

    print(f"Wrote {written_rows} rows for {len(uniprot_ids)} unique UniProt IDs.")


if __name__ == "__main__":
    main()
