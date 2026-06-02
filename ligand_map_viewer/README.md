# Ligand map viewer

Interactive Streamlit app for `specific_ligands.csv`: protein/ligand dropdowns, an AgGrid ligand table with 2D SMILES structures (optional row checkboxes), and 3D PDB views (py3Dmol).

Related data files at repo root (from `notebooks/uniprot_pdb_ligand_mapping.ipynb`):

| File | Description |
|------|-------------|
| `query.csv` | Input protein wishlist (client, protein_name, uniprot_id) |
| `all_ligands.csv` | Full PDB/ligand map for all query proteins |
| `specific_ligands.csv` | Ligands unique to a single UniProt — **used by this app** |

## Local run (MaSIF conda env)

```bash
conda activate MaSIF
streamlit run ligand_map_viewer/app.py --server.port 8501 --server.address 0.0.0.0
```

SSH tunnel: `ssh -L 8501:localhost:8501 <host>` → open http://localhost:8501

## Deploy to Streamlit Community Cloud (free)

### Prerequisites

1. Push this repository to **GitHub** (public repo, or your one free private app).
2. Ensure these files are committed:
   - `specific_ligands.csv` (repo root)
   - `packages.txt` (repo root — one apt package name per line, no comments)
   - `ligand_map_viewer/app.py`
   - `ligand_map_viewer/pdb_mapping.py`
   - `ligand_map_viewer/environment.yml`

### Steps

1. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with GitHub.
2. Click **Create app** → **From existing repo**.
3. Select your repo and branch.
4. Set **Main file path** to: `ligand_map_viewer/app.py`
5. **App URL** (optional): e.g. `ligand-map-viewer`
6. Under **Advanced settings** (recommended):
   - **Python version**: 3.11 (or 3.10 if RDKit build fails on 3.11)
7. Click **Deploy**.

Cloud installs dependencies from **`ligand_map_viewer/environment.yml`** (conda-forge RDKit with Cairo) plus **`packages.txt`** at the repo root (X11/Cairo system libs). Local MaSIF installs can still use `requirements.txt` via pip.

Pushes to the deployed branch redeploy automatically.

### Notes

- **Data**: CSV is loaded from the repo root (`specific_ligands.csv`). No secrets required.
- **Network**: The app calls RCSB (`data.rcsb.org`, `files.rcsb.org`) when you select a row.
- **Cache**: PDB files are cached under `data/pdb_cache/` on the container (ephemeral; cleared on redeploy).
- **Cold start**: Apps sleep when idle; first load after idle may take ~30s.
- **Limits**: Free tier has fair-use CPU/memory limits; avoid tight `st.rerun()` loops.

### Troubleshooting deploy

| Issue | Fix |
|-------|-----|
| `ImportError` on `rdkit.Chem.Draw` / `rdMolDraw2D` | Commit root `packages.txt` + `ligand_map_viewer/environment.yml`; reboot app |
| `AttributeError: _ARRAY_API not found` (NumPy 2.x) | Ensure `environment.yml` pins `numpy<2` (pip must not upgrade NumPy) |
| `ModuleNotFoundError: rdkit` | Retry deploy; ensure `environment.yml` is used (not pip-only `requirements.txt`) |
| CSV not found | Commit `specific_ligands.csv` at repo root |
| Build timeout | Pin lighter versions in `requirements.txt` |
| 3D view empty | Check Streamlit logs; RCSB may be rate-limiting |
