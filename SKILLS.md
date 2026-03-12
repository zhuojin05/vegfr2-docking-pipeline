# SKILLS.md — vegfr2-docking-pipeline

Claude Code built-in skills available in this project. Invoke any skill with
its slash command (e.g. `/rdkit`) at the start of a prompt. The skill expands
into a full context-aware assistant for that tool.

---

## `/rdkit` — RDKit cheminformatics

**What it does:** SMILES parsing, descriptor computation, 3D conformer
generation, SDF/MOL I/O, substructure search, fingerprinting.

**When to use it in this pipeline:**
- Stage 1 (`prepare_ligands.py`): SMILES → 3D conformer → SDF
- Computing Lipinski Ro5 properties (MW, logP, HBD, HBA) and Veber rules
  (rotatable bonds, TPSA) in `ligand_properties.csv`
- Validating SMILES before passing to Meeko

**Example:**
```
/rdkit Generate a 3D conformer for sorafenib SMILES and compute MW, logP, TPSA, and rotatable bond count.
```

---

## `/pdb-database` — RCSB PDB fetching

**What it does:** Downloads PDB/mmCIF structures by accession, searches for
structures by UniProt ID or keyword, retrieves metadata and biological assemblies.

**When to use it in this pipeline:**
- Fetching the receptor `4ASD` (VEGFR2, DFG-out) and reference complex `3WZE`
  (VEGFR2 + sorafenib) if `data/raw/` is empty
- Searching for alternative VEGFR2 structures (different inhibitor classes,
  DFG-in vs DFG-out conformations) to expand the virtual screening context
- Checking which HETATM residue name sorafenib uses in 3WZE before running
  `define_binding_site.py`

**Example:**
```
/pdb-database Fetch 3WZE and list all HETATM residue names in the structure.
```

---

## `/pymol` — PyMOL visualisation

**What it does:** Opens PyMOL sessions, generates publication-quality figures,
aligns structures, measures distances, colours by property.

**When to use it in this pipeline:**
- Sanity-checking receptor prep: confirm waters and heteroatoms removed,
  missing loops added by pdbfixer
- Visualising docked poses from `data/results/` overlaid on the 4ASD receptor
- Comparing the top sorafenib docking pose against the 3WZE co-crystal pose
  (RMSD validation)
- Generating figures highlighting key binding residues (Cys919, Glu885,
  Asp1046, Phe1047)

**Example:**
```
/pymol Load VEGFR2.pdbqt and the top sorafenib pose, align to 3WZE, and measure RMSD.
```

---

## `/chembl-database` — ChEMBL bioactivity data

**What it does:** Queries the ChEMBL REST API for compounds, targets, and
bioactivity data (IC50, Ki, Kd); retrieves SMILES for actives.

**When to use it in this pipeline:**
- Finding additional VEGFR2 inhibitors to add to `ligands.csv` beyond the
  three reference compounds
- Retrieving experimental IC50 values for sorafenib, pazopanib, and imatinib
  against VEGFR2 to validate Gnina CNNaffinity rankings
- Enriching `data/results/docking_summary.csv` with experimental data for
  correlation plots

**Example:**
```
/chembl-database Find all ChEMBL compounds with IC50 < 100 nM against VEGFR2 (target CHEMBL279) and return their SMILES.
```

---

## `/biopython` — BioPython structural biology

**What it does:** PDB file parsing with the Bio.PDB module, sequence handling,
pairwise alignment, NCBI Entrez queries, structure superposition.

**When to use it in this pipeline:**
- Parsing 3WZE to extract sorafenib heavy atom coordinates for binding site
  centroid calculation (Stage 3, `define_binding_site.py`)
- Sequence-level checks after pdbfixer repair: confirm no residues are missing
  in the ATP-binding loop region
- NCBI queries for VEGFR2 (KDR) sequence or annotation data

**Example:**
```
/biopython Parse 3WZE.pdb, extract all HETATM atoms for residue SFN, and compute the centroid coordinates.
```

---

## `/alphafold-database` — AlphaFold structure database

**What it does:** Fetches AlphaFold predicted structures from the EBI database
by UniProt accession; retrieves pLDDT confidence scores.

**When to use it in this pipeline:**
- If docking against a VEGFR2 isoform or point mutant not present in the PDB,
  use the AlphaFold model (UniProt P35968) as receptor input
- Comparing AlphaFold model quality (pLDDT) in the ATP-binding pocket against
  the experimental 4ASD structure before deciding which to use
- Academic context: discussing AF2 confidence in the kinase domain vs the
  flexible extracellular region

**Example:**
```
/alphafold-database Fetch the AlphaFold model for human VEGFR2 (UniProt P35968) and report the mean pLDDT for residues 800–1000 (kinase domain).
```

---

## `/uniprot-database` — UniProt protein annotations

**What it does:** Queries UniProt for protein entries, retrieves sequences,
annotated active sites, PTMs, isoforms, and cross-references.

**When to use it in this pipeline:**
- Confirming the residue numbering of key VEGFR2 ATP-pocket residues (Cys919,
  Glu885, Asp1046, Phe1047) relative to the canonical P35968 sequence
- Mapping PDB chain residue numbers to UniProt positions if they differ
  (important for ProLIF interaction labelling in Stage 5)
- Fetching isoform sequences if screening against a splice variant

**Example:**
```
/uniprot-database Look up UniProt P35968 (VEGFR2) and return the annotated active site and binding site residues with their sequence positions.
```

---

## Quick reference

| Slash command | Primary pipeline stage | Key use case |
|---|---|---|
| `/rdkit` | Stage 1 — ligand prep | SMILES → 3D, Ro5/Veber properties |
| `/pdb-database` | Setup / Stage 3 | Fetch 4ASD, 3WZE; find sorafenib HETATM name |
| `/pymol` | Stage 5 — analysis | Visualise poses, RMSD vs co-crystal |
| `/chembl-database` | Input curation / Stage 5 | Expand ligand set, validate rankings vs IC50 |
| `/biopython` | Stage 3 — binding site | Parse PDB, extract ligand centroid |
| `/alphafold-database` | Receptor prep (mutants) | AF2 model if PDB structure unavailable |
| `/uniprot-database` | Stage 5 — interaction labels | Residue numbering, active site annotation |
