"""Streamlit viewer for specific_ligands.csv with 2D/3D structure display."""

from __future__ import annotations

import base64
import io
from pathlib import Path

import pandas as pd
import streamlit as st
from rdkit import Chem
from st_aggrid import AgGrid, GridOptionsBuilder, JsCode

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
ROW_HEIGHT = 240
STRUCTURE_COL_WIDTH = 440

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
    _style_version: int = 2,
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


def _build_aggrid(df: pd.DataFrame):
    display_cols = [c for c in df.columns if c != "structure"] + ["structure"]
    df_show = df[display_cols].copy()
    df_show["structure"] = df_show["structure"].fillna("")

    gb = GridOptionsBuilder.from_dataframe(df_show)
    gb.configure_default_column(filterable=True, sortable=True, resizable=True)
    gb.configure_selection(selection_mode="single", use_checkbox=False)
    gb.configure_grid_options(domLayout="normal", rowHeight=ROW_HEIGHT)

    image_renderer = JsCode(
        """
        class ImgRenderer {
          init(params) {
            this.eGui = document.createElement('div');
            this.eGui.style.display = 'flex';
            this.eGui.style.alignItems = 'center';
            this.eGui.style.height = '100%%';
            if (params.value) {
              const img = document.createElement('img');
              img.src = 'data:image/png;base64,' + params.value;
              img.style.height = '220px';
              img.style.maxWidth = '400px';
              img.style.objectFit = 'contain';
              this.eGui.appendChild(img);
            }
          }
          getGui() { return this.eGui; }
          refresh() { return false; }
        }
        """
    )
    gb.configure_column(
        "structure",
        headerName="2D structure",
        cellRenderer=image_renderer,
        width=STRUCTURE_COL_WIDTH,
    )
    if "smiles" in df_show.columns:
        gb.configure_column("smiles", width=220)

    grid_options = gb.build()
    return AgGrid(
        df_show,
        gridOptions=grid_options,
        update_on="selectionChanged",
        height=min(1200, 80 + len(df_show) * ROW_HEIGHT),
        theme="streamlit",
        allow_unsafe_jscode=True,
        fit_columns_on_grid_load=False,
    )


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

    with col_table:
        st.subheader("Ligand table")
        st.caption(f"{len(df)} rows · {DEFAULT_CSV.name}")
        grid_response = _build_aggrid(df)
        selected = grid_response.get("selected_rows")

    with col_viewer:
        st.subheader("3D structure")
        row = None
        if isinstance(selected, pd.DataFrame) and not selected.empty:
            row = selected.iloc[0]
        elif isinstance(selected, list) and len(selected) > 0:
            row = selected[0]

        if row is None:
            st.info("Select a row in the table to load the PDB structure.")
            return

        if isinstance(row, pd.Series):
            row = row.to_dict()
        pdb_id = row["pdb_id"]
        uniprot_id = row["uniprot_id"]
        ligand_code = row["ligand_code"]

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
