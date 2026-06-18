#!/usr/bin/env python3
"""Query RCSB PDB by (pdb_id, ligand_code) pairs and stream one row per pair to CSV.

Accepts either:
  --input  CSV with columns pdb_id,ligand_code  (batch mode)
  --pdb_id / --ligand_code                       (single-pair mode)

Output columns (same as uniprot_to_pdb_ligand.py):
  uniprot_id, pdb_id, protein_chain, ligand_chain, ligand_code, ligand_name,
  smiles, percent_intracellular

uniprot_id is derived from the PDB polymer entity that owns the chosen protein
chain.  percent_intracellular is computed via DeepTMHMM + EMBOSS needle exactly
as in uniprot_to_pdb_ligand.py.

API responses are cached under data/pdb_api_cache (shared with
uniprot_to_pdb_ligand.py) so re-runs avoid repeated HTTP calls.
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
    "gene_name",
    "recommendedName",
    "scientificName",
    "pdb_id",
    "protein_chain",
    "ligand_chain",
    "ligand_code",
    "ligand_name",
    "smiles",
    "percent_intracellular",
]

DEFAULT_DEEPTMHMM_DIR = Path("/work/upthomae/Meng/hp_list_DeepTMHMM/out/3line")

CacheSource = Literal["pdb"]

_PDB_CACHE: Optional["ApiCache"] = None
_UNIPROT_CACHE: Optional["ApiCache"] = None
_WRITE_LOCKS: Dict[str, threading.Lock] = {}
_WRITE_LOCKS_GUARD = threading.Lock()
_DEEPTMHMM_DIR: Optional[Path] = None
_DEEPTMHMM_ENABLED: bool = True
_3LINE_CACHE: Dict[str, Optional[Tuple[str, str]]] = {}
_3LINE_CACHE_LOCK = threading.Lock()


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

class ApiCache:
    """Filesystem cache for a single API base URL's JSON responses."""

    def __init__(
        self,
        cache_dir: Path,
        base_url: str,
        *,
        enabled: bool = True,
        refresh: bool = False,
    ) -> None:
        self.cache_dir = cache_dir
        self.base_url = base_url.rstrip("/")
        self.enabled = enabled
        self.refresh = refresh
        if enabled:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_json(self, url: str) -> Tuple[Optional[Dict[str, Any]], int]:
        if not self.enabled:
            return self._fetch_live(url)
        cache_path = self._cache_path(url)
        if not self.refresh and cache_path.is_file():
            return self._read_envelope(cache_path)
        body, status = self._fetch_live(url)
        self._write_envelope(cache_path, url, status, body)
        return body, status

    def _cache_path(self, url: str) -> Path:
        rel = url.removeprefix(self.base_url + "/")
        if not rel.endswith(".json"):
            rel = f"{rel}.json"
        return self.cache_dir / rel

    def _read_envelope(self, path: Path) -> Tuple[Optional[Dict[str, Any]], int]:
        with path.open("r", encoding="utf-8") as fh:
            envelope = json.load(fh)
        status = int(envelope.get("status_code", 0))
        body = envelope.get("body")
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


def _init_api_caches(
    pdb_dir: Path,
    uniprot_dir: Path,
    *,
    enabled: bool,
    refresh: bool,
) -> None:
    global _PDB_CACHE, _UNIPROT_CACHE
    _PDB_CACHE = ApiCache(pdb_dir, RCSB_CORE_API, enabled=enabled, refresh=refresh)
    _UNIPROT_CACHE = ApiCache(uniprot_dir, UNIPROT_API, enabled=enabled, refresh=refresh)


def _init_deeptmhmm(deeptmhmm_dir: Path, *, enabled: bool) -> None:
    global _DEEPTMHMM_DIR, _DEEPTMHMM_ENABLED
    _DEEPTMHMM_DIR = deeptmhmm_dir
    _DEEPTMHMM_ENABLED = enabled


def _fetch_pdb_core_json(path_suffix: str) -> Tuple[Optional[Dict[str, Any]], int]:
    url = f"{RCSB_CORE_API}/{path_suffix}"
    if _PDB_CACHE is None:
        return ApiCache._fetch_live(url)
    return _PDB_CACHE.get_json(url)


def _fetch_uniprot_json(uniprot_id: str) -> Tuple[Optional[Dict[str, Any]], int]:
    url = f"{UNIPROT_API}/uniprotkb/{uniprot_id}.json"
    if _UNIPROT_CACHE is None:
        return ApiCache._fetch_live(url)
    return _UNIPROT_CACHE.get_json(url)


def _get_uniprot_metadata(
    uniprot_id: str,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Return (gene_name, recommendedName, scientificName) for a UniProt accession.

    Falls back to submissionNames for unreviewed entries that lack a
    recommendedName (e.g. TrEMBL records).
    """
    data, status = _fetch_uniprot_json(uniprot_id)
    if status != 200 or data is None:
        return None, None, None

    gene_name: Optional[str] = None
    genes = data.get("genes") or []
    if genes:
        gene_name = (genes[0].get("geneName") or {}).get("value")

    recommended_name: Optional[str] = None
    protein_desc = data.get("proteinDescription") or {}

    rec = protein_desc.get("recommendedName") or {}
    full_name = rec.get("fullName") or {}
    recommended_name = full_name.get("value") or None

    if recommended_name is None:
        # Unreviewed entries use submissionNames instead
        submission_names = protein_desc.get("submissionNames") or []
        if submission_names:
            recommended_name = (
                (submission_names[0].get("fullName") or {}).get("value") or None
            )

    scientific_name: Optional[str] = None
    organism = data.get("organism") or {}
    scientific_name = organism.get("scientificName") or None

    return gene_name, recommended_name, scientific_name


# ---------------------------------------------------------------------------
# Polymer map (same as uniprot_to_pdb_ligand.py)
# ---------------------------------------------------------------------------

def _build_polymer_uniprot_maps(
    pdb_id: str, polymer_entity_ids: Sequence[str]
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, str]]:
    """Map polymer entity/asym -> UniProt IDs and label_asym_id -> auth_asym_id."""
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


# ---------------------------------------------------------------------------
# DeepTMHMM + EMBOSS needle (identical to uniprot_to_pdb_ligand.py)
# ---------------------------------------------------------------------------

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
        parsed = (sequence, topology) if len(sequence) == len(topology) else None

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
                "-asequence", str(path_a),
                "-bsequence", str(path_b),
                "-gapopen", "10",
                "-gapextend", "0.5",
                "-outfile", str(path_out),
                "-aformat", "pair",
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
        auth_ids = [str(v) for v in (container.get("auth_asym_ids") or [])]
        asym_ids = [str(v) for v in (container.get("asym_ids") or [])]
        strand_raw = (data.get("entity_poly") or {}).get("pdbx_strand_id") or ""
        strand_ids = [p.strip() for p in str(strand_raw).split(",") if p.strip()]
        if chain not in auth_ids and chain not in asym_ids and chain not in strand_ids:
            continue
        seq = (data.get("entity_poly") or {}).get("pdbx_seq_one_letter_code")
        if seq:
            return "".join(str(seq).split()).upper()
    return None


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

    pdb_sequence = _polymer_sequence_for_auth_chain(pdb_id, protein_chain, polymer_entity_ids)
    if not pdb_sequence:
        return None

    aligned_pdb, aligned_uniprot = _run_needle_alignment(pdb_sequence, uniprot_sequence)
    if aligned_pdb is None or aligned_uniprot is None:
        return None

    return _percent_intracellular_from_alignment(aligned_pdb, aligned_uniprot, topology)


# ---------------------------------------------------------------------------
# SMILES lookup
# ---------------------------------------------------------------------------

def _get_chemcomp_smiles(comp_id: str) -> Optional[str]:
    data, status = _fetch_pdb_core_json(f"chemcomp/{comp_id}")
    if status != 200 or data is None:
        return None
    descriptor = data.get("rcsb_chem_comp_descriptor") or {}
    return descriptor.get("SMILES_stereo") or descriptor.get("SMILES")


# ---------------------------------------------------------------------------
# Instance-scoring helpers
# ---------------------------------------------------------------------------

def _auth_chain_from_instance(instance_data: Dict[str, Any], label_asym_id: str) -> str:
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
    target_asym = neighbor.get("target_asym_id")
    if not target_asym:
        return None
    seq_id = neighbor.get("target_seq_id")
    if seq_id is None:
        seq_id = neighbor.get("target_auth_seq_id")
    if seq_id is None:
        return None
    return target_asym, seq_id


def _score_instance_polymer_neighbors(
    neighbors: Sequence[Dict[str, Any]],
    polymer_asyms: Set[str],
) -> Optional[Tuple[int, str]]:
    """Score one ligand instance by contacts with any polymer chain.

    Returns (total unique polymer residues contacted, best protein label_asym_id)
    or None if no polymer neighbors exist.
    """
    matched_residues: Set[Tuple[str, Any]] = set()
    per_chain: Dict[str, Set[Tuple[str, Any]]] = {}
    chain_order: List[str] = []

    for neighbor in neighbors:
        target_asym = neighbor.get("target_asym_id")
        if target_asym not in polymer_asyms:
            continue
        residue = _residue_key(neighbor)
        if residue is None:
            continue
        matched_residues.add(residue)
        if target_asym not in per_chain:
            per_chain[target_asym] = set()
            chain_order.append(target_asym)
        per_chain[target_asym].add(residue)

    if not matched_residues:
        return None

    best_protein_asym = chain_order[0]
    best_count = len(per_chain[best_protein_asym])
    for asym in chain_order[1:]:
        count = len(per_chain[asym])
        if count > best_count:
            best_count = count
            best_protein_asym = asym

    return len(matched_residues), best_protein_asym


# ---------------------------------------------------------------------------
# Core query: (pdb_id, ligand_code) -> one row dict
# ---------------------------------------------------------------------------

def get_ligand_from_pdb(
    pdb_id: str, ligand_code: str
) -> Optional[Dict[str, Any]]:
    """Fetch ligand data for a specific ligand_code in a PDB entry.

    Scores all instances by polymer-contact count, picks the best, then:
    - derives uniprot_id from the winning protein chain's entity mapping
    - computes percent_intracellular via DeepTMHMM + EMBOSS needle

    Returns a dict or None if the ligand is absent / has no polymer contacts.
    """
    pdb_id = pdb_id.upper()
    ligand_code = ligand_code.upper()

    entry_data, status = _fetch_pdb_core_json(f"entry/{pdb_id}")
    if status != 200 or entry_data is None:
        return None

    try:
        container = entry_data.get("rcsb_entry_container_identifiers") or {}
    except (TypeError, KeyError):
        return None

    nonpolymer_entity_ids = container.get("non_polymer_entity_ids") or []
    polymer_entity_ids = container.get("polymer_entity_ids") or []
    if not nonpolymer_entity_ids:
        return None

    # Build full polymer maps so we can get UniProt IDs for the chosen chain
    _entity_to_uniprot, asym_to_uniprot, label_asym_to_auth = _build_polymer_uniprot_maps(
        pdb_id, polymer_entity_ids
    )
    polymer_asyms: Set[str] = set(asym_to_uniprot.keys())

    for entity_id in nonpolymer_entity_ids:
        entity_data, status = _fetch_pdb_core_json(f"nonpolymer_entity/{pdb_id}/{entity_id}")
        if status != 200 or entity_data is None:
            continue
        try:
            comp_id = (entity_data.get("pdbx_entity_nonpoly") or {}).get("comp_id")
            comp_name = (entity_data.get("pdbx_entity_nonpoly") or {}).get("name")
        except (TypeError, KeyError):
            continue
        if not comp_id or comp_id.upper() != ligand_code:
            continue

        instance_asym_ids = (
            (entity_data.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get(
                "asym_ids"
            )
            or []
        )

        best_instance: Optional[Tuple[int, int, str, str]] = None
        best_instance_data: Optional[Dict[str, Any]] = None
        # (total_residues, -instance_index, protein_label_asym, ligand_label_asym)

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

            score = _score_instance_polymer_neighbors(neighbors, polymer_asyms)
            if score is None:
                continue

            total_residues, protein_label_asym = score
            candidate = (total_residues, -instance_index, protein_label_asym, ligand_asym)
            if best_instance is None or candidate > best_instance:
                best_instance = candidate
                best_instance_data = instance_data

        if best_instance is None or best_instance_data is None:
            return None

        _, _, protein_label_asym, ligand_label_asym = best_instance
        protein_chain = label_asym_to_auth.get(protein_label_asym, protein_label_asym)
        ligand_chain = _auth_chain_from_instance(best_instance_data, ligand_label_asym)

        # Derive uniprot_id from the chosen protein chain
        uniprot_ids = (
            asym_to_uniprot.get(protein_label_asym)
            or asym_to_uniprot.get(protein_chain)
            or []
        )
        uniprot_id: Optional[str] = uniprot_ids[0] if uniprot_ids else None

        # Fetch gene_name, recommendedName, and scientificName from UniProt
        gene_name: Optional[str] = None
        recommended_name: Optional[str] = None
        scientific_name: Optional[str] = None
        if uniprot_id:
            gene_name, recommended_name, scientific_name = _get_uniprot_metadata(uniprot_id)

        # Compute percent_intracellular using the same logic as uniprot_to_pdb_ligand.py
        pct_intra: Optional[float] = None
        if uniprot_id:
            pct_intra = compute_percent_intracellular(
                uniprot_id, pdb_id, protein_chain, polymer_entity_ids
            )

        return {
            "uniprot_id": uniprot_id,
            "gene_name": gene_name,
            "recommendedName": recommended_name,
            "scientificName": scientific_name,
            "ligand_code": comp_id,
            "ligand_name": comp_name,
            "smiles": _get_chemcomp_smiles(comp_id),
            "protein_chain": protein_chain,
            "ligand_chain": ligand_chain,
            "percent_intracellular": pct_intra,
        }

    return None


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def build_row_for_pair(pdb_id: str, ligand_code: str) -> Dict[str, Any]:
    """Return one output row for the given (pdb_id, ligand_code) pair."""
    try:
        result = get_ligand_from_pdb(pdb_id, ligand_code)
    except (TypeError, KeyError, AttributeError) as exc:
        print(f"{pdb_id} {ligand_code}: skipped ({exc})")
        result = None

    if result is None:
        return {
            "uniprot_id": None,
            "gene_name": None,
            "recommendedName": None,
            "scientificName": None,
            "pdb_id": pdb_id,
            "protein_chain": None,
            "ligand_chain": None,
            "ligand_code": ligand_code,
            "ligand_name": None,
            "smiles": None,
            "percent_intracellular": None,
        }

    return {
        "uniprot_id": result.get("uniprot_id"),
        "gene_name": result.get("gene_name"),
        "recommendedName": result.get("recommendedName"),
        "scientificName": result.get("scientificName"),
        "pdb_id": pdb_id,
        "protein_chain": result.get("protein_chain"),
        "ligand_chain": result.get("ligand_chain"),
        "ligand_code": result.get("ligand_code") or ligand_code,
        "ligand_name": result.get("ligand_name"),
        "smiles": result.get("smiles"),
        "percent_intracellular": result.get("percent_intracellular"),
    }


def iter_parallel_results(
    pairs: Sequence[Tuple[str, str]], n_jobs: int
) -> Iterable[Dict[str, Any]]:
    parallel = Parallel(n_jobs=n_jobs, prefer="threads", return_as="generator_unordered")
    generator = parallel(
        delayed(build_row_for_pair)(pdb_id, ligand_code)
        for pdb_id, ligand_code in pairs
    )
    yield from generator


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_pairs_from_csv(path: Path) -> List[Tuple[str, str]]:
    """Read (pdb_id, ligand_code) pairs from a CSV file (deduplicates)."""
    pairs: List[Tuple[str, str]] = []
    seen: Set[Tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            return pairs
        pdb_col = "pdb_id" if "pdb_id" in reader.fieldnames else "pdb"
        for row in reader:
            pdb_id = row.get(pdb_col, "").strip().upper()
            ligand_code = row.get("ligand_code", "").strip().upper()
            if not pdb_id or not ligand_code:
                continue
            key = (pdb_id, ligand_code)
            if key in seen:
                continue
            seen.add(key)
            pairs.append(key)
    return pairs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Query RCSB PDB by (pdb_id, ligand_code) pairs, one row per pair."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input", type=Path,
        help="CSV with pdb_id,ligand_code columns (batch mode).",
    )
    input_group.add_argument(
        "--pdb_id",
        help="Single PDB ID (use together with --ligand_code).",
    )

    parser.add_argument(
        "--ligand_code",
        help="Ligand component code (required when --pdb_id is used).",
    )
    parser.add_argument("--output", required=True, type=Path, help="Output CSV path.")
    parser.add_argument(
        "--n-jobs", type=int, default=-1,
        help="joblib worker count (default: -1 uses all cores).",
    )
    parser.add_argument(
        "--pdb-cache-dir", type=Path, default=Path("data/pdb_api_cache"),
        help="Directory for cached RCSB PDB Core JSON (default: data/pdb_api_cache).",
    )
    parser.add_argument(
        "--uniprot-cache-dir", type=Path, default=Path("data/uniprot_api_cache"),
        help="Directory for cached UniProt JSON (default: data/uniprot_api_cache).",
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Do not read or write API cache.",
    )
    parser.add_argument(
        "--refresh-cache", action="store_true",
        help="Ignore existing cache files and refetch from network.",
    )
    parser.add_argument(
        "--deeptmhmm-dir", type=Path, default=DEFAULT_DEEPTMHMM_DIR,
        help=f"Directory with DeepTMHMM .3line files (default: {DEFAULT_DEEPTMHMM_DIR}).",
    )
    parser.add_argument(
        "--no-deeptmhmm", action="store_true",
        help="Skip percent_intracellular computation (column will be empty).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.pdb_id is not None and not args.ligand_code:
        raise SystemExit("--ligand_code is required when --pdb_id is used")

    _init_api_caches(
        args.pdb_cache_dir,
        args.uniprot_cache_dir,
        enabled=not args.no_cache,
        refresh=args.refresh_cache,
    )
    _init_deeptmhmm(args.deeptmhmm_dir, enabled=not args.no_deeptmhmm)

    if args.input is not None:
        pairs = read_pairs_from_csv(args.input)
        print(f"Input CSV:  {args.input}")
    else:
        pairs = [(args.pdb_id.strip().upper(), args.ligand_code.strip().upper())]
        print(f"Single pair: {pairs[0]}")

    print(f"Output CSV: {args.output}")
    if args.no_cache:
        print("API cache: disabled")
    else:
        print(f"PDB cache:     {args.pdb_cache_dir}")
        print(f"UniProt cache: {args.uniprot_cache_dir}")
        if args.refresh_cache:
            print("API cache: refresh (overwrite on fetch)")
    if args.no_deeptmhmm:
        print("DeepTMHMM: disabled")
    else:
        print(f"DeepTMHMM: {args.deeptmhmm_dir}")
    print(f"Pairs to resolve: {len(pairs)}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    written = 0

    with args.output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for row in tqdm(
            iter_parallel_results(pairs, args.n_jobs),
            total=len(pairs),
            desc="Pairs",
            unit="pair",
        ):
            writer.writerow(row)
            fh.flush()
            written += 1

    print(f"Wrote {written} rows for {len(pairs)} pairs.")


if __name__ == "__main__":
    main()
