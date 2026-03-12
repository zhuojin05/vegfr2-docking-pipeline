# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

<!-- project title: vegfr2-docking-pipeline -->

## About this project

Automated, fully scriptable molecular docking pipeline for virtual screening
of small molecule inhibitors against VEGFR2 (PDB: 4ASD), a validated
anticancer target. Replaces a GUI-dependent AutoDock Tools + VINA undergraduate
practical (Imperial College London, Protein Science course) with a reproducible
Python pipeline using Meeko, Gnina, and ProLIF.

Reference co-crystal structure: 3WZE (VEGFR2 + sorafenib), used for automatic
binding site extraction and pose validation.

---

## Developer profile

- BSc Biochemistry, Year 2, Imperial College London
- Comfortable with Python; background in computational biology and protein design
- Familiar with: RDKit basics, PDB file format, conda environments, git
- Working on related projects: BAGEL Monte Carlo protein design framework
  (protein-nucleotide interactions), DINOv2 breast ultrasound classification
- Prefers code with explanatory comments on non-obvious biology/cheminformatics
  choices

---

## Hardware and environment

- macOS Apple Silicon (arm64) — this affects every dependency decision
- CPU only, no CUDA GPU
- Managed with uv (`uv sync --extra dev`) — no conda needed
- Working directory: `~/Documents/GitHub/vegfr2-docking-pipeline`
- Python version: 3.10 (pinned in `.python-version`, uv downloads it automatically)

**Always check the uv venv is active before debugging any import or dependency
errors.** `which python` should point to `.venv/bin/python`, or prefix commands
with `uv run`.

---

## Project structure

```
vegfr2-docking-pipeline/
├── CLAUDE.md                    # this file
├── README.md
├── LICENSE
├── pyproject.toml               # all Python deps; uv reads this
├── .python-version              # pins Python 3.10 for uv
├── config.yaml                  # all runtime parameters
├── pipeline.py                  # CLI orchestrator (argparse)
├── data/
│   ├── raw/                     # VEGFR2.pdb (4ASD), 3WZE.pdb — do not modify
│   ├── prepared/                # pdbqt files, cleaned PDBs — generated
│   ├── ligands/                 # per-ligand .sdf, .pdb, .pdbqt — generated
│   │   └── ligands.csv          # input: name, smiles columns
│   └── results/                 # docking outputs, analysis, plots — generated
├── pipeline/
│   ├── __init__.py
│   ├── prepare_ligands.py       # Stage 1: SMILES → 3D → pdbqt
│   ├── prepare_receptor.py      # Stage 2: PDB → pdbfixer → pdbqt
│   ├── define_binding_site.py   # Stage 3: auto box from 3WZE ligand
│   ├── run_docking.py           # Stage 4: Gnina + VINA backends
│   └── analyse_results.py       # Stage 5: ProLIF fingerprints + plots
├── tests/
│   ├── test_ligand_prep.py
│   ├── test_receptor_prep.py
│   └── test_binding_site.py
├── notebooks/
│   └── results_exploration.ipynb
└── logs/                        # timestamped run logs — generated
```

**Never modify files in `data/raw/`.** All outputs are written to
`data/prepared/`, `data/ligands/`, or `data/results/`.

---

## Pipeline stages overview

| Stage | Script | Input | Output |
|-------|--------|-------|--------|
| 1 | `prepare_ligands.py` | `ligands.csv` (SMILES) | `.sdf`, `.pdb`, `.pdbqt` per ligand + `ligand_properties.csv` |
| 2 | `prepare_receptor.py` | `data/raw/VEGFR2.pdb` | `data/prepared/VEGFR2_clean.pdb`, `VEGFR2.pdbqt` |
| 3 | `define_binding_site.py` | `data/raw/3WZE.pdb` | `data/prepared/binding_site.json`, `config_template.txt` |
| 4 | `run_docking.py` | pdbqt files + binding_site.json | poses `.pdbqt` + scores CSV per ligand |
| 5 | `analyse_results.py` | poses + receptor | ProLIF fingerprints, plots, `docking_summary.csv` |

Run all stages sequentially with:
```bash
uv run python pipeline.py run-all --input data/ligands/ligands.csv
```

---

## Key dependencies and Apple Silicon notes

### Install via uv sync (all Python packages)
As of mid-2025, arm64 PyPI wheels exist for all required packages:
- `rdkit` — arm64 wheels on PyPI
- `pdbfixer` / `openmm` — arm64 wheels from OpenMM project (v8.x+)
- `meeko`, `prolif` — pure Python, no issues
- `vina` — no arm64 wheel on PyPI; source build is also broken (setup.py hardcodes
  `-std=c++11`, incompatible with Homebrew Boost 1.86+). A patched arm64 wheel
  (compiled with `-std=c++14`, Boost 1.90) is included in `wheels/` and referenced
  via `[tool.uv.sources]` in `pyproject.toml`. `uv sync` uses it automatically.
  If you ever need to rebuild it:
    ```bash
    brew install boost swig
    # patch /tmp/vina-1.2.7/setup.py: -std=c++11 → -std=c++14
    CONDA_DEFAULT_ENV=base CONDA_PREFIX=/opt/homebrew uv build /tmp/vina-1.2.7 \
        --wheel --out-dir wheels/
    ```
- `gemmi` — required by meeko>=0.5 (added to explicit deps in pyproject.toml)
- `biopython` — used in Stage 3 (`Bio.PDB.PDBParser`) to parse 3WZE and extract ligand coordinates

Run `uv sync --extra dev` to install everything.

### Gnina on Apple Silicon
Gnina requires CUDA 12.0+, which is unavailable on macOS. There is no native
arm64 binary and no Homebrew formula (`brew install gnina` does not exist;
`apt-get` is Linux-only). Rosetta 2 alone is insufficient — the CUDA runtime
is still unavailable on macOS regardless of emulation layer.

The only working option is Docker:

1. Install Docker Desktop for Mac (arm64 native):
   https://www.docker.com/products/docker-desktop/

2. Pull the gnina image (x86_64; runs via Rosetta 2 on Apple Silicon):
   ```bash
   docker pull gnina/gnina:latest
   ```

3. Set backend in `config.yaml`:
   ```yaml
   docking:
     backend: gnina_docker
   ```

**Why Rosetta 2 + Docker works here:**
gnina's `--no_gpu` flag disables CUDA initialisation entirely. Docker Desktop
emulates x86_64 Linux via Rosetta 2, so gnina runs as a CPU-only process.
Performance: ~2–4× slower than native Linux; accurate CNN scoring is preserved.
Expect ~10–30 minutes per ligand on Apple Silicon CPU (x86_64 emulation + CNN
scoring). For faster local runs use `backend: vina` (default).

**Always flag Apple Silicon compatibility issues explicitly. Never silently
work around them — document the workaround in a comment and in the README.**

### Checking environment health
```bash
uv run python -c "import rdkit; print('rdkit ok')"
uv run python -c "import meeko; print('meeko ok')"
uv run python -c "import prolif; print('prolif ok')"
uv run python -c "from vina import Vina; print('vina ok')"
uv run python -c "import pdbfixer; print('pdbfixer ok')"
which gnina && gnina --version
```

---

## Docking backends

### Gnina (preferred)
- CNN-based scoring trained on PDBbind — more accurate than VINA empirical
- CPU mode is slow: use `--num_modes 10` and `--exhaustiveness 8`
- Output per pose: CNNscore, CNNaffinity, minimizedAffinity
- Poses saved as `.pdbqt`

### AutoDock VINA (Python bindings, comparison backend)
- Empirical scoring function, faster on CPU than Gnina
- Used for score comparison against Gnina
- Same receptor/ligand/box inputs as Gnina
- Switch via `config.yaml`: `docking.backend: vina`

---

## Binding site

Extracted automatically from the sorafenib ligand in 3WZE:
- Sorafenib HETATM residue name in 3WZE: **`BAX`** (chain A, residue 1201) — confirmed by `grep HETATM data/raw/3WZE.pdb`. Code falls back to `SFN/SOR/SF4/SRF` if BAX is not found.
- Box centre: centroid of all heavy atoms of the BAX residue
- Box size: (max_coords − min_coords) + 2 × padding (default 10 Å each side)
- Parameters stored in: `data/prepared/binding_site.json`

**Do not use the ICL practical centroid (1.323, 6.779, 9.145) as a sanity check.**
Those values are in the **4ASD coordinate frame**, not 3WZE's — the two structures
are independently deposited and have different origins. The extracted centroid from
3WZE will differ and is correct.

---

## config.yaml reference

```yaml
docking:
  backend: gnina          # gnina or vina
  num_modes: 10
  exhaustiveness: 8
  padding: 10.0           # Angstroms added to each side of binding site box

paths:
  raw_data: data/raw
  prepared: data/prepared
  ligands: data/ligands
  results: data/results
  logs: logs

receptor:
  pdb: VEGFR2.pdb
  reference_complex: 3WZE.pdb
  ph: 7.4
  remove_heteroatoms: true
  remove_waters: true
```

---

## Test ligands

Located at `data/ligands/ligands.csv`:

| name | smiles |
|------|--------|
| sorafenib | `C(F)(F)(F)OC1=CC=C(C=C1)NC(=O)C2=CC=C(O2)CNC(=O)C3=CC=CN=C3Cl` |
| pazopanib | `CC1=C(C=C(C=C1)NC2=NC=CC(=N2)N(C)C3=CC4=NN(C(=C4C=C3)C)C)S(=O)(=O)N` |
| imatinib | `CC1=CC=C(C=C1)NC(=O)C2=CC=C(C=N2)CNC(=O)C3=CC(=CC=C3)CN4CCN(CC4)C` |

Sorafenib is the reference compound (present in 3WZE) — its top docking pose
should recapitulate the co-crystal binding mode. Use RMSD against 3WZE
sorafenib as a validation check for the pipeline.

---

## Code conventions

- Type hints on all function signatures
- Docstrings on all public functions (NumPy style)
- Per-ligand failures should log a warning and continue — never abort the batch
- All file paths resolved relative to project root using `pathlib.Path`
- Timestamps in log filenames: `logs/pipeline_YYYYMMDD_HHMMSS.log`
- No hardcoded absolute paths anywhere in pipeline/ scripts
- Config always loaded from `config.yaml` via the orchestrator and passed down;
  scripts do not load config independently

---

## Git conventions

- Commit after each working stage
- Commit message format: `Stage N: brief description`
- Never commit: `data/prepared/`, `data/results/`, `logs/`, `*.pdbqt`,
  `__pycache__/`, `.DS_Store`
- Always commit: `data/ligands/ligands.csv`, `config.yaml`, `pyproject.toml`,
  `.python-version`, all `pipeline/` scripts, `CLAUDE.md`

---

## Testing

Tests use pytest and must not require Gnina or VINA to be installed.
Test ligand: aspirin (`CC(=O)Oc1ccccc1C(=O)O`)
Run with:
```bash
uv run pytest tests/ -v                              # all tests
uv run pytest tests/test_ligand_prep.py -v           # single file
uv run pytest tests/test_ligand_prep.py::test_name   # single test
```

Tests cover Stages 1–3 only (ligand prep, receptor prep, binding site
extraction). Stage 4–5 tests are integration tests requiring docking binaries
and are skipped in CI.

---

## Common failure modes and fixes

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `ImportError: No module named 'rdkit'` | uv venv not active | Use `uv run python` or activate `.venv` |
| `pdbfixer` import fails | uv sync not run | `uv sync` |
| Meeko fails on ligand | Rotatable bond detection issue | Check SMILES validity in RDKit first |
| `gnina not found in PATH` on macOS | No native macOS binary; CUDA required | Use `backend: gnina_docker` (Docker) or `backend: vina` |
| 3WZE sorafenib not found | Wrong HETATM residue name | `grep HETATM data/raw/3WZE.pdb` to find actual name |
| VINA box too small | Padding insufficient | Increase `docking.padding` in config.yaml |
| ProLIF no interactions found | Pose PDB missing hydrogens | Ensure hydrogens added before ProLIF analysis |
| `mk_prepare_receptor.py` exits 0 but no `.pdbqt` written | Missing `-p` flag (Meeko ≥0.5 requires it to enable PDBQT output) | Ensure `-p` is passed; Meeko's default output is JSON/GPF only |
| `vina` build fails: "Boost library location was not found!" | vina setup.py never reads BOOST_ROOT; also incompatible with Boost ≥1.86 (c++11 vs c++14). Pre-built wheel in `wheels/` bypasses this. | Use the wheel in `wheels/` via `[tool.uv.sources]`; see CLAUDE.md for rebuild instructions |

---

## Scientific context (for generating comments and documentation)

- VEGFR2 (Vascular Endothelial Growth Factor Receptor 2): receptor tyrosine
  kinase driving tumour angiogenesis; approved inhibitors include sorafenib,
  pazopanib, sunitinib
- Sorafenib binds the ATP pocket in the DFG-out (inactive) kinase conformation
- Key binding residues in VEGFR2 ATP pocket: Cys919 (hinge H-bond), Glu885,
  Asp1046 (DFG), Phe1047; ProLIF output should show interactions with these
- CNNaffinity (Gnina) is in predicted pKd units; minimizedAffinity and VINA
  scores are in kcal/mol (more negative = better)
- Lipinski Ro5: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10
- Veber rules: rotatable bonds ≤ 10, TPSA ≤ 140 Å²

---

## Useful references

- Gnina: https://github.com/gnina/gnina
- Meeko: https://github.com/forlilab/Meeko
- ProLIF: https://prolif.readthedocs.io
- AutoDock VINA Python: https://autodock-vina.readthedocs.io
- pdbfixer: https://github.com/openmm/pdbfixer
- 3WZE structure: https://www.rcsb.org/structure/3WZE
- 4ASD structure: https://www.rcsb.org/structure/4ASD