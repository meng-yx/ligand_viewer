#!/bin/bash
# Run the PDB ligand map Dash viewer (requires conda env MaSIF).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

source ~/.bashrc
conda activate MaSIF

exec python apps/pdb_ligand_viewer/app.py "$@"
