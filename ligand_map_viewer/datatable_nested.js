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
    let html =
        '<table class="ligand-child-table display compact stripe hover" style="width:100%">';
    html += "<thead><tr>";
    for (let i = 0; i < headers.length; i++) {
        html += "<th>" + headers[i] + "</th>";
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
        html += "<td>" + esc(lig.label) + "</td>";
        html += "<td>" + esc(lig.pdb_id) + "</td>";
        html += "<td>" + esc(lig.ligand_code) + "</td>";
        html += "<td>" + esc(lig.ligand_name) + "</td>";
        html += "<td>" + esc(lig.formula) + "</td>";
        html += "<td>" + esc(lig.mw) + "</td>";
        html += "<td>" + esc(lig.num_carbon) + "</td>";
        html += "<td>" + structureCell + "</td>";
        html += "</tr>";
    }
    html += "</tbody></table>";
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
    if (window.__ligandNestedTableWired === tableId) {
        return;
    }
    window.__ligandNestedTableWired = tableId;

    const $table = window.jQuery("#" + tableId);

    $table.on("click", "tbody td.dt-control", function (e) {
        e.stopPropagation();
        const tr = window.jQuery(this).closest("tr");
        const row = table.row(tr);
        if (row.child.isShown()) {
            row.child.hide();
            tr.removeClass("dt-hasChild");
        } else {
            const html = window.commonJsFunctions.formatLigandChild(row.data());
            row.child(html).show();
            tr.addClass("dt-hasChild");
        }
    });

    $table.on("click", ".ligand-select-row", function (e) {
        e.stopPropagation();
        const el = window.jQuery(this);
        const ligand = {
            uniprot_id: el.data("uniprotId"),
            pdb_id: el.data("pdbId"),
            ligand_code: el.data("ligandCode"),
            label: el.data("label"),
            ligand_name: el.data("ligandName"),
            formula: el.data("formula"),
            mw: el.data("mw"),
            num_carbon: el.data("numCarbon"),
        };
        if (!ligand.pdb_id) {
            return;
        }
        window.publicJsFunctions.setComponentValue({
            data: ligand,
            eventType: "ligandClick",
        });
        $table.find(".ligand-select-row").removeClass("ligand-selected");
        el.addClass("ligand-selected");
    });
}

function renderDtControl(data, type, row, meta) {
    if (type === "display" && meta.row === 0) {
        window.commonJsFunctions.wireNestedLigandTable();
    }
    return "";
}
