"""RCSB PDB API helpers and structure filtering for ligand_map_viewer."""

from __future__ import annotations

import io
from pathlib import Path
from typing import Iterable

import requests
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa

RCSB_CORE_API = "https://data.rcsb.org/rest/v1/core"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

PROTEIN_COLORS = [
    "cyanCarbon",
    "greenCarbon",
    "magentaCarbon",
    "yellowCarbon",
    "whiteCarbon",
]

LIGAND_COLOR = "greyCarbon"


def _get_entry_container(pdb_id: str) -> dict | None:
    try:
        response = requests.get(f"{RCSB_CORE_API}/entry/{pdb_id.upper()}", timeout=30)
        if response.status_code != 200:
            return None
        return response.json().get("rcsb_entry_container_identifiers") or {}
    except (requests.RequestException, ValueError, TypeError, KeyError):
        return None


def _build_polymer_uniprot_maps(pdb_id: str, polymer_entity_ids: Iterable) -> tuple[dict, dict]:
    """Map polymer entity_id and asym_id -> UniProt accessions for a PDB entry."""
    entity_to_uniprot: dict[str, list] = {}
    asym_to_uniprot: dict[str, list] = {}

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

        for asym_id in (container.get("asym_ids") or []) + (container.get("auth_asym_ids") or []):
            asym_to_uniprot[asym_id] = uniprot_ids

    return entity_to_uniprot, asym_to_uniprot


def _neighbor_uniprot_ids(neighbor: dict, entity_to_uniprot: dict, asym_to_uniprot: dict) -> set[str]:
    """Resolve a target neighbor record to UniProt accessions."""
    uniprot_ids: set[str] = set()

    target_entity_id = neighbor.get("target_entity_id")
    if target_entity_id is not None:
        for uid in entity_to_uniprot.get(str(target_entity_id), []):
            uniprot_ids.add(uid.upper())

    target_asym_id = neighbor.get("target_asym_id")
    if target_asym_id is not None:
        for uid in asym_to_uniprot.get(target_asym_id, []):
            uniprot_ids.add(uid.upper())

    return uniprot_ids


def _chains_for_uniprot(asym_to_uniprot: dict, uniprot_id: str) -> set[str]:
    chains: set[str] = set()
    for asym_id, uids in asym_to_uniprot.items():
        if uniprot_id.upper() in {u.upper() for u in uids}:
            chains.add(asym_id)
    return chains


def _pdb_chain_id(pdb_id: str, asym_id: str) -> str:
    """Map RCSB asym_id to auth_asym_id used in downloadable PDB files."""
    try:
        response = requests.get(
            f"{RCSB_CORE_API}/nonpolymer_entity_instance/{pdb_id}/{asym_id}",
            timeout=30,
        )
        if response.status_code == 200:
            container = (
                response.json().get("rcsb_nonpolymer_entity_instance_container_identifiers")
                or {}
            )
            auth = container.get("auth_asym_id")
            if auth:
                return auth
    except (requests.RequestException, ValueError, TypeError, KeyError):
        pass
    return asym_id


def resolve_chains_and_ligands(
    pdb_id: str,
    uniprot_id: str,
    ligand_code: str,
) -> tuple[list[str], list[str], list[str]]:
    """
    Resolve protein and ligand chain IDs for a PDB row using RCSB neighbor annotations.

    Returns:
        protein_chain_ids, ligand_chain_ids, warnings
    """
    pdb_id = pdb_id.upper()
    uniprot_id = uniprot_id.upper()
    ligand_code = ligand_code.upper()
    warnings: list[str] = []

    container = _get_entry_container(pdb_id)
    if not container:
        raise ValueError(f"Could not fetch RCSB entry metadata for {pdb_id}")

    polymer_entity_ids = container.get("polymer_entity_ids") or []
    nonpolymer_entity_ids = container.get("non_polymer_entity_ids") or []
    entity_to_uniprot, asym_to_uniprot = _build_polymer_uniprot_maps(pdb_id, polymer_entity_ids)

    protein_chains: set[str] = set()
    ligand_chains: set[str] = set()
    used_neighbor_logic = False

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
        if not comp_id or comp_id.upper() != ligand_code:
            continue

        instance_asym_ids = (
            (entity_data.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("asym_ids")
            or []
        )

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

            contacting_protein: set[str] = set()
            for neighbor in neighbors:
                if uniprot_id not in _neighbor_uniprot_ids(neighbor, entity_to_uniprot, asym_to_uniprot):
                    continue
                target_asym = neighbor.get("target_asym_id")
                if target_asym:
                    contacting_protein.add(target_asym)

            if contacting_protein:
                used_neighbor_logic = True
                ligand_chains.add(_pdb_chain_id(pdb_id, asym_id))
                protein_chains.update(
                    {_pdb_chain_id(pdb_id, c) if c else c for c in contacting_protein}
                )

    if not protein_chains:
        protein_chains = _chains_for_uniprot(asym_to_uniprot, uniprot_id)
        if protein_chains:
            warnings.append(
                "No neighbor-linked protein chain; using all UniProt-matching chains."
            )

    if not ligand_chains:
        if used_neighbor_logic:
            warnings.append("Ligand instances found but none contact the query UniProt.")
        else:
            warnings.append(
                "No neighbor-linked ligand instance; will use distance fallback on PDB file."
            )

    if not protein_chains:
        raise ValueError(
            f"No polymer chain maps to UniProt {uniprot_id} in PDB {pdb_id}"
        )

    return sorted(protein_chains), sorted(ligand_chains), warnings


def download_pdb(pdb_id: str, cache_dir: Path) -> Path:
    """Download PDB file to cache_dir and return local path."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.upper()
    path = cache_dir / f"{pdb_id}.pdb"
    if path.exists():
        return path

    url = RCSB_PDB_URL.format(pdb_id=pdb_id)
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    path.write_text(response.text)
    return path


def _expand_chain_id_set(chain_ids: Iterable[str]) -> set[str]:
    return {c.strip() for c in chain_ids if c}


class _ChainSubsetSelect(Select):
    def __init__(self, chain_ids: set[str]):
        self.chain_ids = chain_ids

    def accept_chain(self, chain) -> bool:
        return chain.id in self.chain_ids

    def accept_residue(self, residue) -> bool:
        return True

    def accept_atom(self, atom) -> bool:
        return True


def _min_distance_between_residues(res_a, res_b) -> float:
    min_dist = float("inf")
    for atom_a in res_a:
        for atom_b in res_b:
            d = atom_a - atom_b
            dist = d.norm()
            if dist < min_dist:
                min_dist = dist
    return min_dist


def _ligand_chains_by_distance(
    structure,
    protein_chain_ids: set[str],
    ligand_code: str,
    cutoff: float = 5.0,
) -> set[str]:
    protein_residues = []
    for model in structure:
        for chain in model:
            if chain.id not in protein_chain_ids:
                continue
            for residue in chain:
                if is_aa(residue, standard=True):
                    protein_residues.append(residue)

    if not protein_residues:
        return set()

    ligand_chains: set[str] = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname.strip().upper() != ligand_code:
                    continue
                for prot_res in protein_residues:
                    if _min_distance_between_residues(residue, prot_res) <= cutoff:
                        ligand_chains.add(chain.id)
                        break
    return ligand_chains


def filter_structure_to_pdb_string(
    pdb_path: Path,
    protein_chain_ids: list[str],
    ligand_chain_ids: list[str],
    ligand_code: str,
    warnings: list[str] | None = None,
) -> tuple[str, list[str], list[str], list[str]]:
    """
    Filter a PDB to protein + ligand chains. Returns (pdb_string, protein_ids, ligand_ids, warnings).
    """
    warnings = list(warnings or [])
    ligand_code = ligand_code.upper()
    protein_set = _expand_chain_id_set(protein_chain_ids)
    ligand_set = _expand_chain_id_set(ligand_chain_ids)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("entry", str(pdb_path))

    if not ligand_set:
        ligand_set = _ligand_chains_by_distance(structure, protein_set, ligand_code)
        if ligand_set:
            warnings.append(
                f"Distance fallback: kept ligand chain(s) {sorted(ligand_set)} within 5 Å of protein."
            )
        else:
            for model in structure:
                for chain in model:
                    for residue in chain:
                        if residue.resname.strip().upper() == ligand_code:
                            ligand_set.add(chain.id)
                            break
            if ligand_set:
                warnings.append(
                    f"No contact-based ligand chain; kept all {ligand_code} chains: {sorted(ligand_set)}."
                )

    keep_chains = protein_set | ligand_set
    if not keep_chains:
        raise ValueError("No chains selected for output structure")

    # Map requested chain IDs to IDs present in the file (asym vs auth)
    file_chain_ids = set()
    for model in structure:
        for chain in model:
            file_chain_ids.add(chain.id)

    resolved_keep: set[str] = set()
    for cid in keep_chains:
        if cid in file_chain_ids:
            resolved_keep.add(cid)

    if not resolved_keep:
        raise ValueError(
            f"None of the resolved chain IDs {sorted(keep_chains)} exist in {pdb_path.name}. "
            f"File has chains: {sorted(file_chain_ids)}"
        )

    io_buf = io.StringIO()
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(io_buf, select=_ChainSubsetSelect(resolved_keep))
    pdb_string = io_buf.getvalue()

    resolved_protein = sorted(resolved_keep & protein_set)
    resolved_ligand = sorted(resolved_keep & ligand_set)
    return pdb_string, resolved_protein, resolved_ligand, warnings


def build_py3dmol_html(
    pdb_string: str,
    protein_chain_ids: list[str],
    ligand_chain_ids: list[str],
    ligand_code: str,
    width: int = 700,
    height: int = 600,
) -> str:
    """Build py3Dmol viewer HTML for filtered structure."""
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_string, "pdb")

    # Cartoon/line per chain (protein and any chain that also hosts the ligand).
    all_chains: list[str] = []
    seen: set[str] = set()
    for chain_id in list(protein_chain_ids) + list(ligand_chain_ids):
        if chain_id not in seen:
            seen.add(chain_id)
            all_chains.append(chain_id)

    for idx, chain_id in enumerate(all_chains):
        color = PROTEIN_COLORS[idx % len(PROTEIN_COLORS)]
        view.addStyle(
            {"chain": chain_id},
            {
                "cartoon": {"colorscheme": color},
                "line": {"colorscheme": color},
            },
        )

    # Only the ligand HET — not every residue on the ligand's chain (e.g. shared chain E).
    view.addStyle(
        {"resn": ligand_code.upper()},
        {"stick": {"colorscheme": LIGAND_COLOR}},
    )

    view.zoomTo()
    return view._make_html()  # noqa: SLF001
