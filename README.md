# VEGFR2 Docking Pipeline

Automated, reproducible molecular docking pipeline for virtual screening of small-molecule
inhibitors against **VEGFR2** (PDB: 4ASD), a validated anti-cancer target.

Replaces the GUI-dependent AutoDock Tools + VINA undergraduate practical
(Imperial College London, Protein Science course) with a scriptable Python pipeline
using Meeko, Gnina, and ProLIF.

---

## Scientific background

**VEGFR2** (Vascular Endothelial Growth Factor Receptor 2) is a receptor tyrosine kinase
that drives tumour angiogenesis. Inhibiting VEGFR2 blocks the formation of new blood vessels
that tumours need to grow. Approved VEGFR2 inhibitors include sorafenib, pazopanib, and sunitinib.

The reference structure **3WZE** is a co-crystal of VEGFR2 bound to sorafenib, which occupies
the ATP-binding pocket in the **DFG-out (inactive)** kinase conformation. Key binding residues:
- **Cys919** — hinge hydrogen bond (backbone NH)
- **Glu885** — upper hydrophobic pocket
- **Asp1046 / Phe1047** — DFG motif; sorafenib's urea group contacts Asp1046

Docking validation: the top sorafenib pose should recapitulate the co-crystal binding mode
(RMSD < 2 Å against 3WZE sorafenib; interactions with Cys919 in ProLIF output).

---

## Prerequisites

- **macOS Apple Silicon (arm64)** — see dependency notes below
- [uv](https://docs.astral.sh/uv/) — single tool for Python + package management
- ~4 GB disk space (environment + PDB files)

---

## Quick start

### Step 1 — install uv (if not already installed)

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Step 2 — install all Python dependencies

```bash
# uv reads .python-version (3.10) and downloads Python automatically if needed.
# vina is built from source automatically (configured in pyproject.toml).
uv sync --extra dev
```

### Step 3 — gnina (not a Python package — separate install)

```bash
brew install gnina           # if an arm64 binary is available
# or: see Apple Silicon section below for fallback options
```

### Step 4 — verify

```bash
uv run python -c "import rdkit; print('rdkit ok')"
uv run python -c "import meeko; print('meeko ok')"
uv run python -c "import prolif; print('prolif ok')"
uv run python -c "from vina import Vina; print('vina ok')"
uv run python -c "import pdbfixer; print('pdbfixer ok')"
which gnina && gnina --version || echo "gnina not found — see below"
```

### Step 5 — run the pipeline

```bash
uv run python pipeline.py run-all --input data/ligands/ligands.csv
```

Or run stages individually:

```bash
uv run python pipeline.py prepare-ligands --input data/ligands/ligands.csv
uv run python pipeline.py prepare-receptor
uv run python pipeline.py define-site
uv run python pipeline.py dock --all          # or: --ligand sorafenib
uv run python pipeline.py analyse
```

---

## Gnina on Apple Silicon (arm64)

Gnina is **not a Python package** (not on PyPI) — it must be installed as a binary
separately from the Python environment. Resolution order:

1. **brew install gnina** (if an arm64 build is available):
   ```bash
   brew install gnina
   gnina --version
   ```

2. **Rosetta 2 fallback** (x86_64 emulation, ~2–3× slower):
   ```bash
   arch -x86_64 brew install gnina
   gnina --version
   ```
   If using Rosetta 2, note this in your results — it affects performance benchmarks
   but not docking accuracy.

3. Download a pre-built binary from https://github.com/gnina/gnina/releases.

**Never silently work around this — document the workaround in your results.**

---

## Dependency notes

All Python dependencies are managed by `uv` via `pyproject.toml`. As of mid-2025,
`rdkit`, `pdbfixer`, and `openmm` all have arm64 PyPI wheels, so no conda is needed.

| Package | Install method | arm64 note |
|---|---|---|
| rdkit | `uv sync` (PyPI) | arm64 wheels available |
| pdbfixer | `uv sync` (PyPI) | arm64 wheels available via OpenMM project |
| openmm | `uv sync` (PyPI) | arm64 wheels available since v8.x |
| meeko | `uv sync` (PyPI) | no issues |
| vina | `uv sync` (source build) | `no-binary = ["vina"]` in pyproject.toml handles this |
| prolif | `uv sync` (PyPI) | no issues |
| gnina | `brew install gnina` or manual binary | **not a Python package** |

---

## Pipeline stages

| Stage | Script | Input | Output |
|---|---|---|---|
| 1 | `pipeline/prepare_ligands.py` | `ligands.csv` (SMILES) | `.sdf`, `.pdb`, `.pdbqt` per ligand + `ligand_properties.csv` |
| 2 | `pipeline/prepare_receptor.py` | `data/raw/4ASD.pdb` | `data/prepared/4ASD_clean.pdb`, `4ASD.pdbqt` |
| 3 | `pipeline/define_binding_site.py` | `data/raw/3WZE.pdb` | `data/prepared/binding_site.json`, `config_template.txt` |
| 4 | `pipeline/run_docking.py` | pdbqt files + `binding_site.json` | poses `.pdbqt` + scores CSV per ligand |
| 5 | `pipeline/analyse_results.py` | poses + receptor | ProLIF fingerprints, plots, `docking_summary.csv` |

---

## Example output

After a successful run, `data/results/docking_summary.csv` will contain:

| ligand | gnina_best_score | vina_best_score | top_interacting_residues | lipinski_pass |
|---|---|---|---|---|
| sorafenib | 0.92 | -9.8 | CYS919; GLU885; ASP1046 | True |
| pazopanib | 0.88 | -9.1 | CYS919; PHE1047 | True |
| imatinib | 0.79 | -8.3 | CYS919 | True |

Plots in `data/results/plots/`:
- `score_vs_rank.png` — score decay across pose ranks
- `interaction_heatmap.png` — ProLIF fingerprint heatmap
- `score_comparison.png` — Gnina vs VINA best scores

---

## GUI workflow vs this pipeline

| Feature | GUI workflow | This pipeline |
|---|---|---|
| Reproducibility | Manual, undocumented steps | Fully scripted, config-driven |
| Batch screening | One ligand at a time | Automated batch via CSV |
| Scoring | VINA empirical only | Gnina CNN + VINA comparison |
| Interaction analysis | Manual inspection | ProLIF fingerprints |
| Binding site | Manual box placement | Automatic from co-crystal |
| Logging | None | Timestamped log files |

---

## Running tests

```bash
pytest tests/ -v
```

Tests cover Stages 1–3. Stages 4–5 require docking binaries and are skipped in CI.

---

## Citations

- **Gnina**: McNutt et al., *J. Cheminform.* 2021, 13, 43. https://github.com/gnina/gnina
- **Meeko**: Forli lab. https://github.com/forlilab/Meeko
- **ProLIF**: Bouysset & Fiorucci, *J. Cheminform.* 2021, 13, 72. https://prolif.readthedocs.io
- **AutoDock VINA**: Eberhardt et al., *J. Chem. Inf. Model.* 2021, 61, 3891.
- **RDKit**: Landrum, G. RDKit: Open-source cheminformatics. https://www.rdkit.org
- **PDBFixer / OpenMM**: Eastman et al., *PLOS Comput. Biol.* 2017, 13, e1005659.
- **3WZE**: Gajiwala et al., *PNAS* 2009 (VEGFR2 + sorafenib). https://www.rcsb.org/structure/3WZE
- **4ASD**: https://www.rcsb.org/structure/4ASD
