#!/usr/bin/env python3
"""Build UniProt -> PDB ligand mapping and stream results to CSV.

API responses are cached under data/uniprot_api_cache and data/pdb_api_cache so
reruns (e.g. to add CSV columns from cached JSON) avoid repeated HTTP calls.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Literal, Optional, Sequence, Set, Tuple

import requests
from joblib import Parallel, delayed
from tqdm import tqdm

RCSB_CORE_API = "https://data.rcsb.org/rest/v1/core"
UNIPROT_API = "https://rest.uniprot.org"

OUTPUT_COLUMNS = [
    "uniprot_id",
    "pdb_id",
    "protein_chain",
    "ligand_chain",
    "ligand_code",
    "ligand_name",
    "smiles",
    "percent_intracellular",
]

DEFAULT_DEEPTMHMM_DIR = Path("/work/upthomae/Meng/hp_list_DeepTMHMM/out/3line")

CacheSource = Literal["uniprot", "pdb"]

# Set in main(); used by fetch helpers (safe across joblib threads).
_API_CACHE: Optional["ApiCache"] = None
_WRITE_LOCKS: Dict[str, threading.Lock] = {}
_WRITE_LOCKS_GUARD = threading.Lock()
_DEEPTMHMM_DIR: Optional[Path] = None
_DEEPTMHMM_ENABLED: bool = True
_3LINE_CACHE: Dict[str, Optional[Tuple[str, str]]] = {}
_3LINE_CACHE_LOCK = threading.Lock()


class ApiCache:
    """Filesystem cache for UniProt and RCSB PDB Core API JSON responses."""

    def __init__(
        self,
        uniprot_dir: Path,
        pdb_dir: Path,
        *,
        enabled: bool = True,
        refresh: bool = False,
    ) -> None:
        self.uniprot_dir = uniprot_dir
        self.pdb_dir = pdb_dir
        self.enabled = enabled
        self.refresh = refresh
        if enabled:
            self.uniprot_dir.mkdir(parents=True, exist_ok=True)
            self.pdb_dir.mkdir(parents=True, exist_ok=True)

    def get_json(self, url: str, source: CacheSource) -> Tuple[Optional[Dict[str, Any]], int]:
        """Return (parsed JSON body or None, HTTP status). Uses cache when enabled."""
        if not self.enabled:
            return self._fetch_live(url)

        cache_path = self._cache_path(url, source)
        if not self.refresh and cache_path.is_file():
            return self._read_envelope(cache_path)

        body, status = self._fetch_live(url)
        self._write_envelope(cache_path, url, status, body)
        return body, status

    def _cache_path(self, url: str, source: CacheSource) -> Path:
        if source == "uniprot":
            # https://rest.uniprot.org/uniprotkb/P12345.json -> uniprotkb/P12345.json
            rel = url.removeprefix(f"{UNIPROT_API}/")
            return self.uniprot_dir / rel

        # https://data.rcsb.org/rest/v1/core/entry/8VLB -> entry/8VLB.json
        rel = url.removeprefix(f"{RCSB_CORE_API}/")
        if not rel.endswith(".json"):
            rel = f"{rel}.json"
        return self.pdb_dir / rel

    def _read_envelope(self, path: Path) -> Tuple[Optional[Dict[str, Any]], int]:
        with path.open("r", encoding="utf-8") as handle:
            envelope = json.load(handle)
        status = int(envelope.get("status_code", 0))
        body = envelope.get("body")
        if body is None:
            return None, status
        return body, status

    def _write_envelope(
        self,
        path: Path,
        url: str,
        status_code: int,
        body: Optional[Dict[str, Any]],
    ) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        envelope = {
            "url": url,
            "status_code": status_code,
            "fetched_at": datetime.now(timezone.utc).isoformat(),
            "body": body,
        }
        payload = json.dumps(envelope, separators=(",", ":"))

        lock = _lock_for(path)
        with lock:
            if path.is_file() and not self.refresh:
                return
            tmp = path.with_suffix(path.suffix + ".tmp")
            tmp.write_text(payload, encoding="utf-8")
            tmp.replace(path)

    @staticmethod
    def _fetch_live(url: str) -> Tuple[Optional[Dict[str, Any]], int]:
        try:
            response = requests.get(url, timeout=30)
            status = response.status_code
            if status == 200:
                return response.json(), status
            return None, status
        except (requests.RequestException, ValueError):
            return None, 0


def _lock_for(path: Path) -> threading.Lock:
    key = str(path.resolve())
    with _WRITE_LOCKS_GUARD:
        if key not in _WRITE_LOCKS:
            _WRITE_LOCKS[key] = threading.Lock()
        return _WRITE_LOCKS[key]


def _init_api_cache(
    uniprot_dir: Path,
    pdb_dir: Path,
    *,
    enabled: bool,
    refresh: bool,
) -> None:
    global _API_CACHE
    _API_CACHE = ApiCache(uniprot_dir, pdb_dir, enabled=enabled, refresh=refresh)


def _init_deeptmhmm(deeptmhmm_dir: Path, *, enabled: bool) -> None:
    global _DEEPTMHMM_DIR, _DEEPTMHMM_ENABLED
    _DEEPTMHMM_DIR = deeptmhmm_dir
    _DEEPTMHMM_ENABLED = enabled


def _fetch_uniprot_json(uniprot_id: str) -> Tuple[Optional[Dict[str, Any]], int]:
    url = f"{UNIPROT_API}/uniprotkb/{uniprot_id}.json"
    if _API_CACHE is None:
        return ApiCache._fetch_live(url)
    return _API_CACHE.get_json(url, "uniprot")


def _fetch_pdb_core_json(path_suffix: str) -> Tuple[Optional[Dict[str, Any]], int]:
    url = f"{RCSB_CORE_API}/{path_suffix}"
    if _API_CACHE is None:
        return ApiCache._fetch_live(url)
    return _API_CACHE.get_json(url, "pdb")


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
            code = value.split(",", 1)[0].strip().upper()
            if code:
                codes.add(code)
    return codes


def get_pdb_accessions(uniprot_id: str) -> List[str]:
    """Return unique PDB IDs cross-referenced for a UniProt accession."""
    data, status = _fetch_uniprot_json(uniprot_id)
    if status != 200 or data is None:
        if status and status != 200:
            print(f"{uniprot_id}: UniProt fetch failed (HTTP {status})")
        elif status == 0:
            print(f"{uniprot_id}: UniProt fetch failed (network/parse error)")
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
    data, status = _fetch_pdb_core_json(f"chemcomp/{comp_id}")
    if status != 200 or data is None:
        return None
    descriptor = data.get("rcsb_chem_comp_descriptor") or {}
    return descriptor.get("SMILES_stereo") or descriptor.get("SMILES")


def _build_polymer_uniprot_maps(
    pdb_id: str, polymer_entity_ids: Sequence[str]
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, str]]:
    """Map polymer entity/asym -> UniProt and label_asym_id -> auth_asym_id."""
    entity_to_uniprot: Dict[str, List[str]] = {}
    asym_to_uniprot: Dict[str, List[str]] = {}
    label_asym_to_auth: Dict[str, str] = {}

    for entity_id in polymer_entity_ids:
        data, status = _fetch_pdb_core_json(f"polymer_entity/{pdb_id}/{entity_id}")
        if status != 200 or data is None:
            continue
        try:
            container = data.get("rcsb_polymer_entity_container_identifiers") or {}
        except (TypeError, KeyError):
            continue
        uniprot_ids = container.get("uniprot_ids") or []
        entity_to_uniprot[str(entity_id)] = uniprot_ids

        asym_ids = container.get("asym_ids") or []
        auth_asym_ids = container.get("auth_asym_ids") or []
        if len(asym_ids) == len(auth_asym_ids):
            for label_asym, auth_asym in zip(asym_ids, auth_asym_ids):
                label_asym_to_auth[label_asym] = auth_asym
                asym_to_uniprot[label_asym] = uniprot_ids
                asym_to_uniprot[auth_asym] = uniprot_ids
        else:
            for asym_id in asym_ids + auth_asym_ids:
                label_asym_to_auth.setdefault(asym_id, asym_id)
                asym_to_uniprot[asym_id] = uniprot_ids

    return entity_to_uniprot, asym_to_uniprot, label_asym_to_auth


def _polymer_sequence_for_auth_chain(
    pdb_id: str, protein_chain: str, polymer_entity_ids: Sequence[str]
) -> Optional[str]:
    """Return one-letter polymer sequence for an auth/label/strand chain ID."""
    chain = protein_chain.strip()
    for entity_id in polymer_entity_ids:
        data, status = _fetch_pdb_core_json(f"polymer_entity/{pdb_id}/{entity_id}")
        if status != 200 or data is None:
            continue
        container = data.get("rcsb_polymer_entity_container_identifiers") or {}
        auth_ids = [str(value) for value in (container.get("auth_asym_ids") or [])]
        asym_ids = [str(value) for value in (container.get("asym_ids") or [])]
        strand_raw = (data.get("entity_poly") or {}).get("pdbx_strand_id") or ""
        strand_ids = [part.strip() for part in str(strand_raw).split(",") if part.strip()]
        if chain not in auth_ids and chain not in asym_ids and chain not in strand_ids:
            continue
        seq = (data.get("entity_poly") or {}).get("pdbx_seq_one_letter_code")
        if seq:
            return "".join(str(seq).split()).upper()
    return None


def _load_deeptmhmm_3line(uniprot_id: str) -> Optional[Tuple[str, str]]:
    """Load (uniprot_sequence, topology) from DeepTMHMM .3line file."""
    if not _DEEPTMHMM_ENABLED or _DEEPTMHMM_DIR is None:
        return None

    uid = uniprot_id.upper()
    with _3LINE_CACHE_LOCK:
        if uid in _3LINE_CACHE:
            return _3LINE_CACHE[uid]

    path = _DEEPTMHMM_DIR / f"{uid}.3line"
    if not path.is_file():
        with _3LINE_CACHE_LOCK:
            _3LINE_CACHE[uid] = None
        return None

    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        parsed: Optional[Tuple[str, str]] = None
    else:
        sequence = "".join(lines[1].split()).upper()
        topology = lines[2].strip()
        if len(sequence) != len(topology):
            parsed = None
        else:
            parsed = (sequence, topology)

    with _3LINE_CACHE_LOCK:
        _3LINE_CACHE[uid] = parsed
    return parsed


def _parse_emboss_pair(path: Path) -> Tuple[Optional[str], Optional[str]]:
    """Parse EMBOSS needle pair output into two gapped alignment strings."""
    aligned_a: List[str] = []
    aligned_b: List[str] = []
    block_index = 0
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith(" "):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            fragment = parts[2]
            if block_index % 2 == 0:
                aligned_a.append(fragment)
            else:
                aligned_b.append(fragment)
            block_index += 1

    if not aligned_a or not aligned_b:
        return None, None
    return "".join(aligned_a), "".join(aligned_b)


def _run_needle_alignment(seq_a: str, seq_b: str) -> Tuple[Optional[str], Optional[str]]:
    """Global-align seq_a (PDB chain) to seq_b (UniProt) with EMBOSS needle."""
    if not shutil.which("needle"):
        return None, None

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        path_a = tmp_path / "pdb_chain.fa"
        path_b = tmp_path / "uniprot.fa"
        path_out = tmp_path / "align.pair"
        path_a.write_text(f">pdb_chain\n{seq_a}\n", encoding="utf-8")
        path_b.write_text(f">uniprot\n{seq_b}\n", encoding="utf-8")

        result = subprocess.run(
            [
                "needle",
                "-asequence",
                str(path_a),
                "-bsequence",
                str(path_b),
                "-gapopen",
                "10",
                "-gapextend",
                "0.5",
                "-outfile",
                str(path_out),
                "-aformat",
                "pair",
                "-auto",
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not path_out.is_file():
            return None, None
        return _parse_emboss_pair(path_out)


def _percent_intracellular_from_alignment(
    aligned_pdb: str, aligned_uniprot: str, topology: str
) -> Optional[float]:
    """Fraction of doubly-aligned columns where UniProt topology is intracellular (I)."""
    if len(aligned_pdb) != len(aligned_uniprot):
        return None

    uniprot_index = 0
    aligned_columns = 0
    intracellular_columns = 0
    for pdb_char, uniprot_char in zip(aligned_pdb, aligned_uniprot):
        if pdb_char != "-" and uniprot_char != "-":
            aligned_columns += 1
            if uniprot_index < len(topology) and topology[uniprot_index] == "I":
                intracellular_columns += 1
        if uniprot_char != "-":
            uniprot_index += 1

    if aligned_columns == 0:
        return None
    return intracellular_columns / aligned_columns


def compute_percent_intracellular(
    uniprot_id: str,
    pdb_id: str,
    protein_chain: Optional[str],
    polymer_entity_ids: Sequence[str],
) -> Optional[float]:
    """Map PDB chain residues to DeepTMHMM topology via Needle; return fraction in 'I'."""
    if not protein_chain:
        return None

    deeptmhmm = _load_deeptmhmm_3line(uniprot_id)
    if deeptmhmm is None:
        return None
    uniprot_sequence, topology = deeptmhmm

    pdb_sequence = _polymer_sequence_for_auth_chain(
        pdb_id, protein_chain, polymer_entity_ids
    )
    if not pdb_sequence:
        return None

    aligned_pdb, aligned_uniprot = _run_needle_alignment(pdb_sequence, uniprot_sequence)
    if aligned_pdb is None or aligned_uniprot is None:
        return None

    return _percent_intracellular_from_alignment(aligned_pdb, aligned_uniprot, topology)


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


def _auth_chain_from_instance(
    instance_data: Dict[str, Any], label_asym_id: str
) -> str:
    """Return auth_asym_id for an entity instance, falling back to label asym_id."""
    for container_key in (
        "rcsb_nonpolymer_entity_instance_container_identifiers",
        "rcsb_polymer_entity_instance_container_identifiers",
    ):
        container = instance_data.get(container_key) or {}
        auth = container.get("auth_asym_id")
        if auth:
            return auth
    return label_asym_id


def _residue_key(neighbor: Dict[str, Any]) -> Optional[Tuple[str, Any]]:
    """Unique protein residue identifier from a target-neighbor record."""
    target_asym = neighbor.get("target_asym_id")
    if not target_asym:
        return None
    seq_id = neighbor.get("target_seq_id")
    if seq_id is None:
        seq_id = neighbor.get("target_auth_seq_id")
    if seq_id is None:
        return None
    return target_asym, seq_id


def _score_ligand_instance_neighbors(
    neighbors: Sequence[Dict[str, Any]],
    uniprot_id: str,
    entity_to_uniprot: Dict[str, List[str]],
    asym_to_uniprot: Dict[str, List[str]],
) -> Optional[Tuple[int, str]]:
    """
    Score one ligand instance by query-matching neighbor residues.

    Returns (total unique residues across query-matching chains, best protein label_asym_id)
    or None if no query-matching neighbors.
    """
    matched_residues: Set[Tuple[str, Any]] = set()
    per_chain: Dict[str, Set[Tuple[str, Any]]] = {}
    chain_order: List[str] = []

    for neighbor in neighbors:
        if uniprot_id not in _neighbor_uniprot_ids(
            neighbor, entity_to_uniprot, asym_to_uniprot
        ):
            continue
        residue = _residue_key(neighbor)
        if residue is None:
            continue
        target_asym = residue[0]
        matched_residues.add(residue)
        if target_asym not in per_chain:
            per_chain[target_asym] = set()
            chain_order.append(target_asym)
        per_chain[target_asym].add(residue)

    if not matched_residues:
        return None

    best_protein_asym = chain_order[0]
    best_count = len(per_chain[best_protein_asym])
    for target_asym in chain_order[1:]:
        count = len(per_chain[target_asym])
        if count > best_count:
            best_count = count
            best_protein_asym = target_asym

    return len(matched_residues), best_protein_asym


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

    Returns a list of dicts with: ligand_code, ligand_name, smiles, protein_chain,
    ligand_chain (auth_asym_id with label-asym fallback).
    """
    pdb_id = pdb_id.upper()
    uniprot_id = uniprot_id.upper()

    entry_data, status = _fetch_pdb_core_json(f"entry/{pdb_id}")
    if status != 200 or entry_data is None:
        return []

    try:
        container = entry_data.get("rcsb_entry_container_identifiers") or {}
    except (TypeError, KeyError):
        return []

    nonpolymer_entity_ids = container.get("non_polymer_entity_ids") or []
    polymer_entity_ids = container.get("polymer_entity_ids") or []
    if not nonpolymer_entity_ids:
        return []

    entity_to_uniprot, asym_to_uniprot, label_asym_to_auth = _build_polymer_uniprot_maps(
        pdb_id, polymer_entity_ids
    )

    ligand_rows: List[Dict[str, Optional[str]]] = []
    seen_comp_ids: Set[str] = set()

    for entity_id in nonpolymer_entity_ids:
        entity_data, status = _fetch_pdb_core_json(
            f"nonpolymer_entity/{pdb_id}/{entity_id}"
        )
        if status != 200 or entity_data is None:
            continue
        try:
            comp_id = (entity_data.get("pdbx_entity_nonpoly") or {}).get("comp_id")
            comp_name = (entity_data.get("pdbx_entity_nonpoly") or {}).get("name")
        except (TypeError, KeyError):
            continue
        if not comp_id or comp_id in seen_comp_ids:
            continue

        if filter_by_subject_of_investigation and not _is_subject_of_investigation(
            entity_data
        ):
            continue

        instance_asym_ids = (
            (entity_data.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get(
                "asym_ids"
            )
            or []
        )

        best_instance: Optional[Tuple[int, int, str, str]] = None
        best_instance_data: Optional[Dict[str, Any]] = None
        # (total_residue_count, -instance_index, protein_label_asym, ligand_label_asym)

        for instance_index, ligand_asym in enumerate(instance_asym_ids):
            instance_data, status = _fetch_pdb_core_json(
                f"nonpolymer_entity_instance/{pdb_id}/{ligand_asym}"
            )
            if status != 200 or instance_data is None:
                continue
            try:
                neighbors = instance_data.get("rcsb_target_neighbors") or []
            except (TypeError, KeyError):
                continue

            score = _score_ligand_instance_neighbors(
                neighbors, uniprot_id, entity_to_uniprot, asym_to_uniprot
            )
            if score is None:
                continue

            total_residues, protein_label_asym = score
            candidate = (total_residues, -instance_index, protein_label_asym, ligand_asym)
            if best_instance is None or candidate > best_instance:
                best_instance = candidate
                best_instance_data = instance_data

        if best_instance is None or best_instance_data is None:
            continue

        _, _, protein_label_asym, ligand_label_asym = best_instance
        protein_chain = label_asym_to_auth.get(protein_label_asym, protein_label_asym)
        ligand_chain = _auth_chain_from_instance(best_instance_data, ligand_label_asym)

        seen_comp_ids.add(comp_id)
        ligand_rows.append(
            {
                "ligand_code": comp_id,
                "ligand_name": comp_name,
                "smiles": _get_chemcomp_smiles(comp_id),
                "protein_chain": protein_chain,
                "ligand_chain": ligand_chain,
                "percent_intracellular": compute_percent_intracellular(
                    uniprot_id, pdb_id, protein_chain, polymer_entity_ids
                ),
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
        except (TypeError, KeyError, AttributeError) as exc:
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
                    "protein_chain": None,
                    "ligand_chain": None,
                    "ligand_code": None,
                    "ligand_name": None,
                    "smiles": None,
                    "percent_intracellular": None,
                }
            )
            continue

        for ligand in filtered_ligands:
            output_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "pdb_id": pdb_id,
                    "protein_chain": ligand.get("protein_chain"),
                    "ligand_chain": ligand.get("ligand_chain"),
                    "ligand_code": ligand.get("ligand_code"),
                    "ligand_name": ligand.get("ligand_name"),
                    "smiles": ligand.get("smiles"),
                    "percent_intracellular": ligand.get("percent_intracellular"),
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
    parser.add_argument(
        "--uniprot-cache-dir",
        type=Path,
        default=Path("data/uniprot_api_cache"),
        help="Directory for cached UniProt JSON (default: data/uniprot_api_cache).",
    )
    parser.add_argument(
        "--pdb-cache-dir",
        type=Path,
        default=Path("data/pdb_api_cache"),
        help="Directory for cached RCSB PDB Core JSON (default: data/pdb_api_cache).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Do not read or write API cache (always fetch from network).",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Ignore existing cache files and refetch from network.",
    )
    parser.add_argument(
        "--deeptmhmm-dir",
        type=Path,
        default=DEFAULT_DEEPTMHMM_DIR,
        help=(
            "Directory with DeepTMHMM .3line files "
            f"(default: {DEFAULT_DEEPTMHMM_DIR})."
        ),
    )
    parser.add_argument(
        "--no-deeptmhmm",
        action="store_true",
        help="Skip percent_intracellular computation (column will be empty).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    blacklist_path = Path(args.blacklist) if args.blacklist else None

    _init_api_cache(
        args.uniprot_cache_dir,
        args.pdb_cache_dir,
        enabled=not args.no_cache,
        refresh=args.refresh_cache,
    )
    _init_deeptmhmm(args.deeptmhmm_dir, enabled=not args.no_deeptmhmm)

    print(f"Input path: {input_path}")
    print(f"Output path: {output_path}")
    print(f"Blacklist path: {blacklist_path}")
    if args.no_cache:
        print("API cache: disabled")
    else:
        print(f"UniProt cache: {args.uniprot_cache_dir}")
        print(f"PDB cache: {args.pdb_cache_dir}")
        if args.refresh_cache:
            print("API cache: refresh (overwrite on fetch)")
    if args.no_deeptmhmm:
        print("DeepTMHMM: disabled")
    else:
        print(f"DeepTMHMM dir: {args.deeptmhmm_dir}")

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
