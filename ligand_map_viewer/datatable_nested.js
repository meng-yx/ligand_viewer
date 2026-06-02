function formatLigandChild(row) {
    const ligands = row.ligands || [];
    if (!ligands.length) {
        return '<div class="ligand-child-empty">No ligands for this protein.</div>';
    }

    const headers = [
        "Label",
        "PDB",
        "Ligand code",
        "Ligand name",
        "Formula",
        "MW",
        "C atoms",
        "2D structure",
    ];
    const headerClasses = [
        "ligand-col-label",
        "ligand-col-pdb",
        "ligand-col-code",
        "ligand-col-name",
        "ligand-col-formula",
        "ligand-col-mw",
        "ligand-col-carbon",
        "ligand-col-structure",
    ];
    let html = '<div class="ligand-child-wrap">';
    html += '<table class="ligand-child-table display compact stripe hover">';
    html += "<thead><tr>";
    for (let i = 0; i < headers.length; i++) {
        html += '<th class="' + headerClasses[i] + '">' + headers[i] + "</th>";
    }
    html += "</tr></thead><tbody>";

    for (let i = 0; i < ligands.length; i++) {
        const lig = ligands[i];
        const esc = function (v) {
            if (v === null || v === undefined) {
                return "";
            }
            return String(v)
                .replace(/&/g, "&amp;")
                .replace(/</g, "&lt;")
                .replace(/"/g, "&quot;");
        };
        let structureCell = "";
        if (lig.structure) {
            structureCell =
                '<img src="data:image/png;base64,' +
                lig.structure +
                '" style="height:120px;max-width:200px;object-fit:contain;" alt="2D structure"/>';
        }
        html +=
            '<tr class="ligand-select-row" ' +
            'data-uniprot-id="' + esc(lig.uniprot_id) + '" ' +
            'data-pdb-id="' + esc(lig.pdb_id) + '" ' +
            'data-ligand-code="' + esc(lig.ligand_code) + '" ' +
            'data-label="' + esc(lig.label) + '" ' +
            'data-ligand-name="' + esc(lig.ligand_name) + '" ' +
            'data-formula="' + esc(lig.formula) + '" ' +
            'data-mw="' + esc(lig.mw) + '" ' +
            'data-num-carbon="' + esc(lig.num_carbon) + '">';
        html += '<td class="ligand-col-label">' + esc(lig.label) + "</td>";
        html += '<td class="ligand-col-pdb">' + esc(lig.pdb_id) + "</td>";
        html += '<td class="ligand-col-code">' + esc(lig.ligand_code) + "</td>";
        html += '<td class="ligand-col-name">' + esc(lig.ligand_name) + "</td>";
        html += '<td class="ligand-col-formula">' + esc(lig.formula) + "</td>";
        html += '<td class="ligand-col-mw">' + esc(lig.mw) + "</td>";
        html += '<td class="ligand-col-carbon">' + esc(lig.num_carbon) + "</td>";
        html += '<td class="ligand-col-structure">' + structureCell + "</td>";
        html += "</tr>";
    }
    html += "</tbody></table></div>";
    return html;
}

function wireNestedLigandTable() {
    if (!window.publicJsFunctions || !window.publicJsFunctions.getDataTableApi) {
        return;
    }
    const table = window.publicJsFunctions.getDataTableApi();
    if (!table) {
        return;
    }
    const tableNode = table.table().node();
    const tableId = tableNode ? tableNode.id : null;
    if (!tableId) {
        return;
    }
    const tableEl = table.table().node();
    if (!tableEl) {
        return;
    }

    // Rebind on every rerender safely (same id, fresh DOM/table instance).
    if (window.__ligandNestedClickHandler) {
        tableEl.removeEventListener("click", window.__ligandNestedClickHandler);
    }

    window.__ligandNestedClickHandler = function (e) {
        const controlCell = e.target.closest("td.dt-control");
        if (controlCell && tableEl.contains(controlCell)) {
            e.stopPropagation();
            const tr = controlCell.closest("tr");
            const row = table.row(tr);
            if (row.child.isShown()) {
                row.child.hide();
                tr.classList.remove("dt-hasChild");
            } else {
                const html = window.commonJsFunctions.formatLigandChild(row.data());
                row.child(html).show();
                tr.classList.add("dt-hasChild");
            }
            return;
        }

        const ligandEl = e.target.closest(".ligand-select-row");
        if (!ligandEl || !tableEl.contains(ligandEl)) {
            return;
        }
        e.stopPropagation();
        const ligand = {
            uniprot_id: ligandEl.dataset.uniprotId,
            pdb_id: ligandEl.dataset.pdbId,
            ligand_code: ligandEl.dataset.ligandCode,
            label: ligandEl.dataset.label,
            ligand_name: ligandEl.dataset.ligandName,
            formula: ligandEl.dataset.formula,
            mw: ligandEl.dataset.mw,
            num_carbon: ligandEl.dataset.numCarbon,
        };
        if (!ligand.pdb_id) {
            return;
        }
        window.publicJsFunctions.setComponentValue({
            data: ligand,
            eventType: "ligandClick",
        });
        const selected = tableEl.querySelectorAll(".ligand-select-row.ligand-selected");
        selected.forEach(function (node) {
            node.classList.remove("ligand-selected");
        });
        ligandEl.classList.add("ligand-selected");
    };

    tableEl.addEventListener("click", window.__ligandNestedClickHandler);
}

function renderDtControl(data, type, row, meta) {
    if (type === "display" && meta.row === 0) {
        window.setTimeout(function () {
            if (
                window.commonJsFunctions &&
                typeof window.commonJsFunctions.wireNestedLigandTable === "function"
            ) {
                window.commonJsFunctions.wireNestedLigandTable();
            }
        }, 0);
    }
    return "";
}
