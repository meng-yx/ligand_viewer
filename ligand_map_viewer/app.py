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
DEFAULT_CSV = REPO_ROOT / "all_ligands.csv"

STRUCTURE_IMG_SIZE = (420, 320)
ROW_HEIGHT = 240
DETAIL_STRUCTURE_COL_WIDTH = 280
SUMMARY_ROW_HEIGHT = 36

ALWAYS_DETAIL_COLUMNS = [
    "uniprot_id",
    "gene_name",
    "recommendedName",
    "pdb_id",
    "ligand_code",
    "ligand_name",
    "smiles",
    "structure",
]
DEFAULT_OPTIONAL_DETAIL_COLUMNS = ["qed"]

DEFAULT_FILTER_RANGES = {
    "mw": (100, None),
    "num_carbon": (5, None),
    "num_N_O": (1, None),
    "qed": (0.0, 1.0),
    "uniprot_id_count": (1, 1),
}



def _smiles_to_base64(smiles: str, size: tuple[int, int] = STRUCTURE_IMG_SIZE) -> str | None:
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
def _cached_smiles_to_base64(smiles: str) -> str | None:
    """Per-SMILES PNG cache; only invoked for rows shown in the detail table."""
    return _smiles_to_base64(smiles)


@st.cache_data(show_spinner=False)
def load_table(csv_path: str) -> pd.DataFrame:
    return pd.read_csv(csv_path)


def attach_structure_column(df: pd.DataFrame) -> pd.DataFrame:
    """Add structure PNG column for visible detail rows only."""
    out = df.copy()
    if "smiles" not in out.columns:
        out["structure"] = ""
        return out

    structures: list[str] = []
    for smi in out["smiles"]:
        if pd.isna(smi) or not str(smi).strip():
            structures.append("")
            continue
        smi_str = str(smi).strip()
        b64 = _cached_smiles_to_base64(smi_str)
        structures.append(b64 or "")
    out["structure"] = structures
    return out


def _pdb_view_url(pdb_id: str) -> str:
    return f"https://www.rcsb.org/3d-view/{str(pdb_id).upper()}"


def _is_liganded(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


def _optional_detail_columns(df: pd.DataFrame) -> list[str]:
    excluded = set(ALWAYS_DETAIL_COLUMNS) | {"is_liganded", "is_valid_ligand"}
    return [c for c in df.columns if c not in excluded]


def _present_values(series: pd.Series) -> pd.Series:
    """Non-missing values (not NaN and not empty/whitespace-only)."""
    empty = series.isna() | series.fillna("").astype(str).str.strip().eq("")
    return series[~empty]


def _optional_column_is_numeric(series: pd.Series) -> bool:
    """True when every present value is numeric; missing values are ignored."""
    values = _present_values(series)
    if values.empty:
        return False
    return bool(pd.to_numeric(values, errors="coerce").notna().all())


def _numeric_ligand_filter_columns(df: pd.DataFrame) -> list[str]:
    """Optional columns that are numeric on all non-missing ligand-coded rows."""
    if "ligand_code" not in df.columns:
        return []
    liganded = df[_is_liganded(df["ligand_code"])]
    numeric_cols: list[str] = []
    for col in _optional_detail_columns(df):
        if col not in liganded.columns:
            continue
        if _optional_column_is_numeric(liganded[col]):
            numeric_cols.append(col)
    return numeric_cols


def _ligand_filter_numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    if "ligand_code" not in df.columns or col not in df.columns:
        return pd.Series(dtype=float)
    liganded = df[_is_liganded(df["ligand_code"])]
    return pd.to_numeric(liganded[col], errors="coerce").dropna()


def _is_integer_column(numeric: pd.Series) -> bool:
    if numeric.empty:
        return False
    return bool((numeric % 1).abs().max() < 1e-9)


def _ligand_filter_bounds(df: pd.DataFrame, col: str) -> tuple[float, float, bool] | None:
    numeric = _ligand_filter_numeric_series(df, col)
    if numeric.empty:
        return None
    min_val, max_val = float(numeric.min()), float(numeric.max())
    return min_val, max_val, _is_integer_column(numeric)


def _resolve_default_slider_range(
    col_name: str,
    min_val: float,
    max_val: float,
    use_int: bool,
) -> tuple[float, float]:
    """Map DEFAULT_FILTER_RANGES onto data bounds; None uses min/max from data."""
    default_lo, default_hi = DEFAULT_FILTER_RANGES.get(col_name, (None, None))
    lo = float(min_val if default_lo is None else default_lo)
    hi = float(max_val if default_hi is None else default_hi)
    lo = max(min_val, min(lo, max_val))
    hi = max(min_val, min(hi, max_val))
    if lo > hi:
        lo, hi = hi, lo
    if use_int:
        lo, hi = int(round(lo)), int(round(hi))
    return float(lo), float(hi)


def render_ligand_filters(df: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Render global ligand-qualification sliders; returns per-column [min, max] ranges."""
    filter_cols = _numeric_ligand_filter_columns(df)
    ranges: dict[str, tuple[float, float]] = {}

    with st.expander("Ligand filters", expanded=True):
        st.caption(
            "Valid ligand = non-empty ligand_code and within all numeric ranges below. "
        )
        if not filter_cols:
            st.info("No numeric optional columns available for ligand filtering.")
            return ranges

        cols_per_row = 3
        for row_start in range(0, len(filter_cols), cols_per_row):
            row_cols = st.columns(cols_per_row)
            for col_name, ui_col in zip(filter_cols[row_start : row_start + cols_per_row], row_cols):
                bounds = _ligand_filter_bounds(df, col_name)
                if bounds is None:
                    continue
                min_val, max_val, use_int = bounds
                if use_int:
                    min_val, max_val = int(round(min_val)), int(round(max_val))
                default_lo, default_hi = _resolve_default_slider_range(
                    col_name, float(min_val), float(max_val), use_int
                )
                with ui_col:
                    if min_val == max_val:
                        label = f"{min_val}" if use_int else f"{min_val:g}"
                        st.caption(f"**{col_name}**: fixed at {label}")
                        ranges[col_name] = (float(min_val), float(max_val))
                    elif use_int:
                        selected = st.slider(
                            col_name,
                            min_value=min_val,
                            max_value=max_val,
                            value=(int(default_lo), int(default_hi)),
                            step=1,
                            key=f"ligand_filter_{col_name}",
                        )
                        ranges[col_name] = (float(selected[0]), float(selected[1]))
                    else:
                        selected = st.slider(
                            col_name,
                            min_value=min_val,
                            max_value=max_val,
                            value=(default_lo, default_hi),
                            key=f"ligand_filter_{col_name}",
                        )
                        ranges[col_name] = (float(selected[0]), float(selected[1]))
    return ranges


def _ranges_cache_key(ranges: dict[str, tuple[float, float]]) -> tuple[tuple[str, float, float], ...]:
    return tuple(sorted((col, lo, hi) for col, (lo, hi) in ranges.items()))


def valid_ligand_mask(df: pd.DataFrame, filter_ranges: dict[str, tuple[float, float]]) -> pd.Series:
    """True when row has ligand_code and passes all numeric ligand filter ranges."""
    if "ligand_code" not in df.columns:
        return pd.Series(False, index=df.index)
    mask = _is_liganded(df["ligand_code"])
    for col, (lo, hi) in filter_ranges.items():
        if col not in df.columns:
            mask &= False
            continue
        numeric = pd.to_numeric(df[col], errors="coerce")
        in_range = numeric.ge(lo) & numeric.le(hi)
        mask &= in_range.fillna(False)
    return mask


@st.cache_data(show_spinner=False)
def build_summary_table(df: pd.DataFrame, filter_ranges_key: tuple[tuple[str, float, float], ...]) -> pd.DataFrame:
    filter_ranges = {col: (lo, hi) for col, lo, hi in filter_ranges_key}
    base = df.copy()
    base["qed_num"] = pd.to_numeric(base["qed"], errors="coerce") if "qed" in base.columns else float("nan")
    base["is_valid_ligand"] = valid_ligand_mask(base, filter_ranges)
    valid_ligands = base[base["is_valid_ligand"]]

    agg: dict = {"pdb_entries": ("pdb_id", "nunique")}
    if "gene_name" in base.columns:
        agg["gene_name"] = ("gene_name", "first")
    if "recommendedName" in base.columns:
        agg["recommendedName"] = ("recommendedName", "first")

    summary = base.groupby("uniprot_id", dropna=False).agg(**agg).reset_index()
    liganded_counts = (
        valid_ligands.groupby("uniprot_id", dropna=False)
        .agg(
            liganded_structures=("pdb_id", "nunique"),
            unique_ligands=("ligand_code", "nunique"),
        )
        .reset_index()
    )
    qed_summary = (
        valid_ligands.groupby("uniprot_id", dropna=False)["qed_num"]
        .median()
        .fillna(0)
        .reset_index(name="median_qed")
    )
    summary = summary.merge(liganded_counts, on="uniprot_id", how="left")
    summary = summary.merge(qed_summary, on="uniprot_id", how="left")
    summary["liganded_structures"] = summary["liganded_structures"].fillna(0).astype(int)
    summary["unique_ligands"] = summary["unique_ligands"].fillna(0).astype(int)
    summary["median_qed"] = summary["median_qed"].fillna(0.0)
    summary = summary.sort_values(["liganded_structures", "pdb_entries", "uniprot_id"], ascending=[False, False, True])
    col_order = [
        "uniprot_id",
        "gene_name",
        "recommendedName",
        "pdb_entries",
        "liganded_structures",
        "unique_ligands",
        "median_qed",
    ]
    return summary[[c for c in col_order if c in summary.columns]]


def _build_summary_aggrid(df: pd.DataFrame):
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_default_column(filterable=True, sortable=True, resizable=True)
    gb.configure_selection(selection_mode="single", use_checkbox=False)
    gb.configure_grid_options(domLayout="normal", rowHeight=SUMMARY_ROW_HEIGHT)
    gb.configure_column("uniprot_id", headerName="UniProt ID", width=120)
    gb.configure_column("gene_name", headerName="Gene", width=110)
    gb.configure_column("recommendedName", headerName="Recommended name", width=260)
    gb.configure_column("pdb_entries", headerName="# PDB entries", width=130)
    gb.configure_column("liganded_structures", headerName="# Liganded structures", width=170)
    gb.configure_column("unique_ligands", headerName="# Unique ligands", width=140)
    gb.configure_column("median_qed", headerName="median_qed", width=110)

    return AgGrid(
        df,
        gridOptions=gb.build(),
        update_on=["selectionChanged"],
        height=min(520, 80 + len(df) * SUMMARY_ROW_HEIGHT),
        theme="streamlit",
        fit_columns_on_grid_load=False,
    )


def _build_detail_aggrid(df: pd.DataFrame, display_cols: list[str]):
    cols = [c for c in display_cols if c in df.columns]
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
    gb.configure_column("uniprot_id", width=110)
    gb.configure_column("gene_name", width=100)
    gb.configure_column("recommendedName", width=220)
    gb.configure_column("ligand_code", width=100)
    gb.configure_column("ligand_name", width=160)
    gb.configure_column("smiles", width=180)
    for optional_col in ["formula", "mw", "qed", "num_carbon", "num_N_O", "uniprot_id_count"]:
        if optional_col in df_show.columns:
            gb.configure_column(optional_col, width=100)
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
    ligand_filter_ranges = render_ligand_filters(df)
    filter_ranges_key = _ranges_cache_key(ligand_filter_ranges)
    summary_df = build_summary_table(df, filter_ranges_key)

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
            detail_df = df[df["uniprot_id"].astype(str) == selected_uniprot].copy()

            optional_cols = _optional_detail_columns(detail_df)
            default_optional = [c for c in DEFAULT_OPTIONAL_DETAIL_COLUMNS if c in optional_cols]

            ctrl_mode, ctrl_cols = st.columns([1, 2], gap="medium")
            with ctrl_mode:
                view_mode = st.segmented_control(
                    "Show",
                    options=["Liganded", "All"],
                    default="Liganded",
                    key="detail_view_mode",
                )
            with ctrl_cols:
                selected_optional_cols = st.multiselect(
                    "Optional columns",
                    options=optional_cols,
                    default=default_optional,
                    help="Always-on columns are fixed; choose additional columns to display.",
                    key="detail_optional_columns",
                )

            if view_mode == "Liganded":
                detail_df = detail_df[valid_ligand_mask(detail_df, ligand_filter_ranges)]

            display_cols = ALWAYS_DETAIL_COLUMNS + selected_optional_cols

            st.caption(
                f"{selected_uniprot} · {len(detail_df)} rows · "
                f"mode: {view_mode}"
            )
            if detail_df.empty:
                st.warning("No rows for this UniProt under the selected mode.")
                detail_selected = None
            else:
                with st.spinner("Generating 2D structures for visible rows…"):
                    detail_df = attach_structure_column(detail_df)
                detail_response = _build_detail_aggrid(detail_df, display_cols=display_cols)
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
