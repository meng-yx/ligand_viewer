"""Streamlit viewer for specific_ligands.csv with 2D/3D structure display."""

from __future__ import annotations

import base64
import io
import json
from pathlib import Path

import pandas as pd
import streamlit as st
from rdkit import Chem
from st_aggrid import AgGrid, JsCode

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


def _structure_image_renderer() -> JsCode:
    return JsCode(
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


def _build_master_df(df: pd.DataFrame) -> pd.DataFrame:
    """One master row per uniprot_id; ligands column is JSON for the detail grid."""
    master_rows: list[dict] = []
    # uniprot_id is the groupby index; sort only by display columns.
    sort_cols = ["label", "protein_name"]
    present = [c for c in sort_cols if c in df.columns]
    grouped = df.groupby("uniprot_id", sort=False)
    if present:
        order = (
            df.groupby("uniprot_id", sort=False)[present]
            .first()
            .sort_values(present)
            .index
        )
    else:
        order = sorted(grouped.groups.keys())

    display_cols = [c for c in df.columns if c != "structure"] + ["structure"]
    for uniprot_id in order:
        group = grouped.get_group(uniprot_id)
        detail = group[display_cols].copy()
        detail["structure"] = detail["structure"].fillna("")
        master_rows.append(
            {
                "uniprot_id": uniprot_id,
                "label": group["label"].iloc[0] if "label" in group.columns else "",
                "protein_name": group["protein_name"].iloc[0]
                if "protein_name" in group.columns
                else "",
                "n_ligand_codes": int(group["ligand_code"].nunique()),
                "ligands": json.dumps(detail.to_dict("records")),
            }
        )
    return pd.DataFrame(master_rows)


def _build_master_detail_grid(df: pd.DataFrame):
    """Collapsed summary per UniProt; expand row for ligand table + 2D structures."""
    master_df = _build_master_df(df)
    image_renderer = _structure_image_renderer()

    detail_column_defs = [
        {"field": "label", "width": 70},
        {"field": "pdb_id", "width": 90},
        {"field": "ligand_code", "width": 100},
        {"field": "ligand_name", "minWidth": 180},
        {"field": "smiles", "width": 220},
        {"field": "formula", "width": 120},
        {"field": "mw", "width": 110},
        {"field": "num_carbon", "headerName": "C atoms", "width": 90},
        {
            "field": "structure",
            "headerName": "2D structure",
            "cellRenderer": image_renderer,
            "width": STRUCTURE_COL_WIDTH,
        },
    ]

    grid_options = {
        "masterDetail": True,
        "detailRowAutoHeight": True,
        "rowSelection": "single",
        "suppressRowClickSelection": False,
        "domLayout": "normal",
        "rowHeight": 44,
        "columnDefs": [
            {
                "field": "uniprot_id",
                "cellRenderer": "agGroupCellRenderer",
                "width": 120,
            },
            {"field": "label", "width": 70},
            {"field": "protein_name", "minWidth": 200, "flex": 1},
            {
                "field": "n_ligand_codes",
                "headerName": "Unique ligands",
                "width": 140,
            },
            {"field": "ligands", "hide": True},
        ],
        "defaultColDef": {
            "filterable": True,
            "sortable": True,
            "resizable": True,
        },
        "detailCellRendererParams": {
            "detailGridOptions": {
                "columnDefs": detail_column_defs,
                "defaultColDef": {
                    "filterable": True,
                    "sortable": True,
                    "resizable": True,
                },
                "rowSelection": "single",
                "domLayout": "normal",
                "rowHeight": ROW_HEIGHT,
            },
            "getDetailRowData": JsCode(
                """function (params) {
                    params.successCallback(JSON.parse(params.data.ligands));
                }"""
            ),
        },
    }

    return AgGrid(
        master_df,
        gridOptions=grid_options,
        update_on="selectionChanged",
        height=min(900, 100 + len(master_df) * 48),
        theme="streamlit",
        allow_unsafe_jscode=True,
        enable_enterprise_modules=True,
        fit_columns_on_grid_load=False,
        key="ligand_master_detail_grid",
    )


def _ligand_row_from_selection(selected) -> dict | None:
    """Extract a ligand detail row (has pdb_id) from AgGrid selection."""
    if selected is None:
        return None
    if isinstance(selected, pd.DataFrame):
        if selected.empty:
            return None
        row = selected.iloc[0]
    elif isinstance(selected, list):
        if not selected:
            return None
        row = selected[0]
    else:
        return None

    if isinstance(row, pd.Series):
        row = row.to_dict()
    if not isinstance(row, dict) or "pdb_id" not in row:
        return None
    return row


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

    with col_table:
        st.subheader("Ligand table")
        st.caption(
            f"{n_proteins} proteins · {len(df)} ligand entries · "
            f"expand a UniProt row for ligand details · {DEFAULT_CSV.name}"
        )
        grid_response = _build_master_detail_grid(df)
        row = _ligand_row_from_selection(grid_response.get("selected_rows"))

    with col_viewer:
        st.subheader("3D structure")

        if row is None:
            st.info(
                "Expand a protein row, then select a ligand in the detail table "
                "to load the PDB structure."
            )
            return
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
