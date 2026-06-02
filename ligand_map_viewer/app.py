"""Streamlit viewer for specific_ligands.csv with 2D/3D structure display."""

from __future__ import annotations

import base64
import io
from pathlib import Path

import pandas as pd
import streamlit as st
from rdkit import Chem

from datatable_ui import ligand_row_from_datatable, render_nested_ligand_table
from pdb_mapping import (
    build_py3dmol_html,
    download_pdb,
    filter_structure_to_pdb_string,
    resolve_chains_and_ligands,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV = REPO_ROOT / "specific_ligands.csv"
PDB_CACHE_DIR = REPO_ROOT / "data" / "pdb_cache"

STRUCTURE_IMG_SIZE = (420, 320)


def _smiles_to_base64(smiles: str | None, size: tuple[int, int] = STRUCTURE_IMG_SIZE) -> str | None:
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return None
    try:
        from rdkit.Chem import Draw
    except ImportError as exc:
        raise ImportError(
            "RDKit Draw failed to load (Cairo/X11 libraries missing). "
            "On Streamlit Cloud, use environment.yml + packages.txt from the repo README."
        ) from exc
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return base64.b64encode(buf.getvalue()).decode("ascii")


@st.cache_data(show_spinner=False)
def load_table(csv_path: str, _structure_img_version: int = 2) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    structures = []
    for smi in df["smiles"]:
        structures.append(_smiles_to_base64(smi))
    df = df.copy()
    df["structure"] = structures
    return df


@st.cache_data(show_spinner="Resolving chains via RCSB…")
def _resolve_cached(pdb_id: str, uniprot_id: str, ligand_code: str) -> tuple[list[str], list[str], list[str]]:
    return resolve_chains_and_ligands(pdb_id, uniprot_id, ligand_code)


@st.cache_data(show_spinner="Downloading PDB…")
def _get_pdb_path(pdb_id: str, cache_dir: str) -> str:
    return str(download_pdb(pdb_id, Path(cache_dir)))


@st.cache_data(show_spinner="Building 3D view…")
def _build_structure_view(
    pdb_id: str,
    uniprot_id: str,
    ligand_code: str,
    cache_dir: str,
    _style_version: int = 3,
) -> tuple[str | None, str, list[str]]:
    """Returns (html, caption, warnings) or raises."""
    protein_chains, ligand_chains, warnings = _resolve_cached(pdb_id, uniprot_id, ligand_code)
    pdb_path = Path(_get_pdb_path(pdb_id, cache_dir))
    pdb_string, protein_chains, ligand_chains, warnings = filter_structure_to_pdb_string(
        pdb_path,
        protein_chains,
        ligand_chains,
        ligand_code,
        warnings,
    )
    html = build_py3dmol_html(pdb_string, protein_chains, ligand_chains, ligand_code)
    caption = (
        f"**{pdb_id}** · protein chains: `{', '.join(protein_chains) or '—'}` · "
        f"ligand chains: `{', '.join(ligand_chains) or '—'}` ({ligand_code})"
    )
    return html, caption, warnings


def main():
    st.set_page_config(page_title="Ligand map viewer", layout="wide")
    st.title("Ligand map viewer")

    if not DEFAULT_CSV.exists():
        st.error(f"CSV not found: {DEFAULT_CSV}")
        st.stop()

    if "viewer_html" not in st.session_state:
        st.session_state.viewer_html = None
        st.session_state.viewer_caption = None
        st.session_state.viewer_warnings = []

    df = load_table(str(DEFAULT_CSV))

    col_table, col_viewer = st.columns([2, 1], gap="medium")

    n_proteins = df["uniprot_id"].nunique()
    ligand_row: dict | None = None

    with col_table:
        st.subheader("Ligand table")
        st.caption(
            f"{n_proteins} proteins · {len(df)} ligand entries · "
            f"expand a UniProt row for ligand details · {DEFAULT_CSV.name}"
        )
        grid_response = render_nested_ligand_table(df)
        ligand_row = ligand_row_from_datatable(grid_response)

    with col_viewer:
        st.subheader("3D structure")

        if ligand_row is None:
            st.info(
                "Expand a protein row, then click a ligand in the nested table "
                "(e.g. pdb_id or ligand_code) to load the PDB structure."
            )
            return

        pdb_id = ligand_row["pdb_id"]
        uniprot_id = ligand_row["uniprot_id"]
        ligand_code = ligand_row["ligand_code"]

        with st.spinner(f"Loading {pdb_id} ({ligand_code})…"):
            try:
                html, caption, warnings = _build_structure_view(
                    str(pdb_id),
                    str(uniprot_id),
                    str(ligand_code),
                    str(PDB_CACHE_DIR),
                )
                st.session_state.viewer_html = html
                st.session_state.viewer_caption = caption
                st.session_state.viewer_warnings = warnings
            except Exception as exc:
                st.session_state.viewer_html = None
                st.session_state.viewer_caption = None
                st.session_state.viewer_warnings = []
                st.error(f"Failed to load structure: {exc}")
                return

        st.html(st.session_state.viewer_html, unsafe_allow_javascript=True)
        st.markdown(st.session_state.viewer_caption)
        for w in st.session_state.viewer_warnings:
            st.warning(w)


if __name__ == "__main__":
    main()
