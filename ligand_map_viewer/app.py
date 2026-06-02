"""Streamlit viewer for specific_ligands.csv with 2D/3D structure display."""

from __future__ import annotations

import base64
import io
from pathlib import Path

import pandas as pd
import streamlit as st
from rdkit import Chem
from st_aggrid import AgGrid, GridOptionsBuilder, JsCode

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV = REPO_ROOT / "specific_ligands.csv"

STRUCTURE_IMG_SIZE = (420, 320)
ROW_HEIGHT = 240
DETAIL_STRUCTURE_COL_WIDTH = 280
SUMMARY_ROW_HEIGHT = 36

DETAIL_COLUMNS = [
    "pdb_id",
    "ligand_code",
    "ligand_name",
    "smiles",
    "formula",
    "mw",
    "num_carbon",
    "structure",
]

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


def _pdb_view_url(pdb_id: str) -> str:
    return f"https://www.rcsb.org/3d-view/{str(pdb_id).upper()}"


def _is_liganded(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


@st.cache_data(show_spinner=False)
def build_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    base = df.copy()
    base["is_liganded"] = _is_liganded(base["ligand_code"])
    liganded = base[base["is_liganded"]]

    agg: dict = {"pdb_entries": ("pdb_id", "nunique")}
    if "protein_name" in base.columns:
        agg["protein_name"] = ("protein_name", "first")

    summary = base.groupby("uniprot_id", dropna=False).agg(**agg).reset_index()
    liganded_counts = (
        liganded.groupby("uniprot_id", dropna=False)
        .agg(
            liganded_structures=("pdb_id", "nunique"),
            unique_ligands=("ligand_code", "nunique"),
        )
        .reset_index()
    )
    summary = summary.merge(liganded_counts, on="uniprot_id", how="left")
    summary["liganded_structures"] = summary["liganded_structures"].fillna(0).astype(int)
    summary["unique_ligands"] = summary["unique_ligands"].fillna(0).astype(int)
    summary = summary.sort_values(["liganded_structures", "pdb_entries", "uniprot_id"], ascending=[False, False, True])
    col_order = ["uniprot_id", "protein_name", "pdb_entries", "liganded_structures", "unique_ligands"]
    return summary[[c for c in col_order if c in summary.columns]]


def _build_summary_aggrid(df: pd.DataFrame):
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_default_column(filterable=True, sortable=True, resizable=True)
    gb.configure_selection(selection_mode="single", use_checkbox=False)
    gb.configure_grid_options(domLayout="normal", rowHeight=SUMMARY_ROW_HEIGHT)
    gb.configure_column("uniprot_id", headerName="UniProt ID", width=120)
    gb.configure_column("protein_name", headerName="Protein", width=140)
    gb.configure_column("pdb_entries", headerName="# PDB entries", width=130)
    gb.configure_column("liganded_structures", headerName="# Liganded structures", width=170)
    gb.configure_column("unique_ligands", headerName="# Unique ligands", width=140)

    return AgGrid(
        df,
        gridOptions=gb.build(),
        update_on=["selectionChanged"],
        height=min(520, 80 + len(df) * SUMMARY_ROW_HEIGHT),
        theme="streamlit",
        fit_columns_on_grid_load=False,
    )


def _build_detail_aggrid(df: pd.DataFrame):
    cols = [c for c in DETAIL_COLUMNS if c in df.columns]
    df_show = df[cols].copy()
    if "structure" in df_show.columns:
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
              img.style.maxWidth = '260px';
              img.style.objectFit = 'contain';
              this.eGui.appendChild(img);
            }
          }
          getGui() { return this.eGui; }
          refresh() { return false; }
        }
        """
    )
    gb.configure_column("pdb_id", width=90)
    gb.configure_column("ligand_code", width=100)
    gb.configure_column("ligand_name", width=160)
    gb.configure_column("smiles", width=180)
    gb.configure_column("formula", width=110)
    gb.configure_column("mw", width=80)
    gb.configure_column("num_carbon", width=90)
    gb.configure_column(
        "structure",
        headerName="2D structure",
        cellRenderer=image_renderer,
        width=DETAIL_STRUCTURE_COL_WIDTH,
    )

    grid_options = gb.build()
    return AgGrid(
        df_show,
        gridOptions=grid_options,
        update_on=["selectionChanged"],
        height=min(1200, 80 + len(df_show) * ROW_HEIGHT),
        theme="streamlit",
        allow_unsafe_jscode=True,
        fit_columns_on_grid_load=False,
    )


def _extract_selected_row(selected) -> dict | None:
    """Return first selected row as a dict across st-aggrid return formats."""
    if isinstance(selected, pd.DataFrame):
        if selected.empty:
            return None
        return selected.iloc[0].to_dict()

    if isinstance(selected, list):
        if not selected:
            return None
        first = selected[0]
        if isinstance(first, pd.Series):
            return first.to_dict()
        if isinstance(first, dict):
            return first

    if isinstance(selected, dict):
        return selected

    return None


def main():
    st.set_page_config(page_title="Ligand map viewer", layout="wide")
    st.title("Ligand map viewer")

    if not DEFAULT_CSV.exists():
        st.error(f"CSV not found: {DEFAULT_CSV}")
        st.stop()

    df = load_table(str(DEFAULT_CSV))

    summary_df = build_summary_table(df)

    col_table, col_viewer = st.columns([3, 2], gap="medium")

    with col_table:
        st.subheader("UniProt summary")
        st.caption(f"{summary_df['uniprot_id'].nunique()} UniProt IDs · source: {DEFAULT_CSV.name}")
        summary_response = _build_summary_aggrid(summary_df)
        summary_selected = _extract_selected_row(summary_response.get("selected_rows"))

        st.divider()
        st.subheader("PDB details")

        if summary_selected is None:
            st.info("Select a UniProt row above.")
            detail_selected = None
        else:
            selected_uniprot = str(summary_selected["uniprot_id"])
            view_mode = st.segmented_control(
                "Show",
                options=["Liganded", "All"],
                default="Liganded",
                key="detail_view_mode",
            )
            detail_df = df[df["uniprot_id"].astype(str) == selected_uniprot].copy()
            if view_mode == "Liganded":
                detail_df = detail_df[_is_liganded(detail_df["ligand_code"])]

            st.caption(
                f"{selected_uniprot} · {len(detail_df)} rows · "
                f"mode: {view_mode}"
            )
            if detail_df.empty:
                st.warning("No rows for this UniProt under the selected mode.")
                detail_selected = None
            else:
                detail_response = _build_detail_aggrid(detail_df)
                detail_selected = _extract_selected_row(detail_response.get("selected_rows"))

    with col_viewer:
        st.subheader("3D structure")
        if detail_selected is None:
            st.info("Select a row in the bottom table to load the PDB structure.")
            return

        pdb_id = str(detail_selected["pdb_id"]).upper()
        view_url = _pdb_view_url(pdb_id)
        st.markdown(f"**PDB:** `{pdb_id}`")
        st.link_button("Open RCSB 3D View", view_url, use_container_width=True)
        st.components.v1.iframe(view_url, height=760, scrolling=True)


if __name__ == "__main__":
    main()
