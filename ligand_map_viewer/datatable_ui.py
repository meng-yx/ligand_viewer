"""DataTables.net nested protein/ligand table for the ligand map viewer."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st
from streamlit_datatables_net import st_datatable, stringify_javascript_function

_JS_PATH = Path(__file__).resolve().parent / "datatable_nested.js"


def _load_js_function(name: str) -> str:
    return stringify_javascript_function(str(_JS_PATH), name)


def _build_nested_rows(df: pd.DataFrame) -> list[dict]:
    """Aggregate rows with nested ligand records for DataTables child rows."""
    rows: list[dict] = []
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

    ligand_fields = [
        "label",
        "uniprot_id",
        "pdb_id",
        "ligand_code",
        "ligand_name",
        "formula",
        "mw",
        "num_carbon",
        "structure",
    ]

    for uniprot_id in order:
        group = grouped.get_group(uniprot_id)
        ligands: list[dict] = []
        for _, lig_row in group.iterrows():
            entry = {f: lig_row[f] for f in ligand_fields if f in lig_row.index}
            structure = entry.get("structure")
            entry["structure"] = structure if structure else ""
            ligands.append(entry)

        rows.append(
            {
                "uniprot_id": uniprot_id,
                "label": group["label"].iloc[0] if "label" in group.columns else "",
                "protein_name": group["protein_name"].iloc[0]
                if "protein_name" in group.columns
                else "",
                "n_ligand_codes": int(group["ligand_code"].nunique()),
                "ligands": ligands,
            }
        )
    return rows


def render_nested_ligand_table(df: pd.DataFrame, *, key: str = "ligand_nested_grid"):
    """Render one expandable DataTables.net table (protein → ligands)."""
    nested_rows = _build_nested_rows(df)
    render_dt_control = _load_js_function("renderDtControl")

    columns = [
        {
            "className": "dt-control",
            "orderable": False,
            "data": None,
            "defaultContent": "",
            "render": render_dt_control,
        },
        {"data": "uniprot_id", "title": "UniProt"},
        {"data": "label", "title": "Label"},
        {"data": "protein_name", "title": "Protein"},
        {"data": "n_ligand_codes", "title": "Unique ligands"},
    ]

    options = {
        "columns": columns,
        "data": nested_rows,
        "order": [[1, "asc"]],
        "scrollX": True,
        "scrollY": "520px",
        "scrollCollapse": True,
        "stateSave": False,
        "paging": False,
        "autoWidth": False,
        "layout": {"topStart": "search"},
        "select": False,
    }

    common_js_functions = {
        "formatLigandChild": _load_js_function("formatLigandChild"),
        "wireNestedLigandTable": _load_js_function("wireNestedLigandTable"),
        "renderDtControl": render_dt_control,
    }

    st.markdown(
        """
        <style>
        div.dt-container div.dt-scroll-body { max-height: 520px !important; }
        div.dt-container table.dataTable { width: 100% !important; }
        .ligand-child-wrap {
            max-height: 460px;
            overflow-y: auto;
            overflow-x: auto;
            margin: 4px 0 0 0;
            padding: 0 0 14px 0;
            box-sizing: border-box;
        }
        tr.dt-hasChild > td { overflow: visible !important; padding-bottom: 10px !important; }
        .ligand-child-table {
            margin: 0;
            font-size: 13px;
            table-layout: fixed;
            min-width: 980px;
            width: max-content !important;
        }
        .ligand-child-table th,
        .ligand-child-table td { white-space: nowrap; }
        .ligand-child-table th.ligand-col-name,
        .ligand-child-table td.ligand-col-name {
            min-width: 260px;
            max-width: 420px;
            white-space: normal;
            word-break: break-word;
        }
        .ligand-select-row { cursor: pointer; }
        .ligand-select-row.ligand-selected { background-color: #ffebee !important; }
        td.dt-control { cursor: pointer; }
        </style>
        """,
        unsafe_allow_html=True,
    )

    return st_datatable(
        nested_rows,
        options=options,
        dt_config="default",
        common_js_functions=common_js_functions,
        override_click_response=True,
        on_select="rerun",
        key=key,
    )


def ligand_row_from_datatable(response) -> dict | None:
    """Parse a ligand dict from st_datatable click / custom component value."""
    if not response:
        return None

    data = None
    if isinstance(response, dict):
        data = response.get("data")
        if data is None and isinstance(response.get("customComponentValue"), dict):
            data = response["customComponentValue"].get("data")

    if isinstance(data, dict) and data.get("pdb_id"):
        return data
    return None
