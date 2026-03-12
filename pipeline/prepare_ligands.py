"""
Stage 1: Prepare ligands from SMILES.

Converts SMILES strings to 3D structures, minimises geometry, writes SDF/PDB/PDBQT
files, and computes physicochemical descriptors (Lipinski Ro5, Veber rules).

Meeko's mk_prepare_ligand.py is called via subprocess to generate PDBQT files
suitable for AutoDock VINA and Gnina.
"""

import logging
import subprocess
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

logger = logging.getLogger(__name__)


def _embed_and_minimise(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Generate a 3D conformer using ETKDGv3 and minimise with MMFF94.

    ETKDGv3 is the recommended conformer generation method as of RDKit 2022+.
    It uses knowledge-based torsion angle distributions for better initial geometry.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule with explicit hydrogens added.

    Returns
    -------
    Chem.Mol or None
        Molecule with 3D coordinates, or None if embedding failed.
    """
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # reproducible conformer generation
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        logger.warning("ETKDGv3 embedding failed")
        return None

    # MMFF94 is a standard force field for drug-like molecules
    ff_result = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    if ff_result == 1:
        logger.warning("MMFF94 minimisation did not converge (maxIters=2000)")
    elif ff_result == -1:
        logger.warning("MMFF94 minimisation failed — using unminimised geometry")

    return mol


def _compute_descriptors(mol: Chem.Mol, name: str) -> dict:
    """
    Compute Lipinski Rule of Five and Veber rule descriptors.

    Lipinski Ro5: MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10
    Veber rules: RotBonds ≤ 10, TPSA ≤ 140 Å²

    VEGFR2 inhibitors like sorafenib (MW 465) and pazopanib (MW 474)
    comply with Lipinski but sit near the upper bounds.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule (may or may not have 3D coords).
    name : str
        Ligand name for logging.

    Returns
    -------
    dict
        Descriptor values and rule compliance flags.
    """
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = Descriptors.TPSA(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    lipinski = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
    veber = (rot_bonds <= 10) and (tpsa <= 140)

    return {
        "name": name,
        "MW": round(mw, 2),
        "MolLogP": round(logp, 3),
        "NumHDonors": hbd,
        "NumHAcceptors": hba,
        "TPSA": round(tpsa, 2),
        "NumRotatableBonds": rot_bonds,
        "Lipinski_pass": lipinski,
        "Veber_pass": veber,
    }


def _run_meeko(sdf_path: Path, pdbqt_path: Path) -> bool:
    """
    Call Meeko's mk_prepare_ligand.py to generate PDBQT from SDF.

    Meeko correctly handles rotatable bonds, atom typing, and partial charges
    for AutoDock VINA / Gnina input. It is preferred over the legacy
    AutoDock Tools prepare_ligand4.py.

    Parameters
    ----------
    sdf_path : Path
        Input SDF file path.
    pdbqt_path : Path
        Output PDBQT file path.

    Returns
    -------
    bool
        True if successful, False otherwise.
    """
    cmd = [
        "mk_prepare_ligand.py",
        "-i", str(sdf_path),
        "-o", str(pdbqt_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.warning(
            "mk_prepare_ligand.py failed for %s: %s",
            sdf_path.stem, result.stderr.strip()
        )
        return False
    return True


def prepare_ligands(csv_path: Path, config: dict) -> pd.DataFrame:
    """
    Stage 1: Convert SMILES to 3D structures and PDBQT files.

    Per-ligand failures log a warning and continue — the batch is never aborted.
    This matches the pipeline design philosophy: partial results are useful,
    and a single bad SMILES should not block screening of other compounds.

    Parameters
    ----------
    csv_path : Path
        CSV with columns 'name' and 'smiles'.
    config : dict
        Pipeline config loaded from config.yaml.

    Returns
    -------
    pd.DataFrame
        Descriptor table for successfully processed ligands.
        Also written to data/ligands/ligand_properties.csv.
    """
    ligands_dir = Path(config["paths"]["ligands"])
    ligands_dir.mkdir(parents=True, exist_ok=True)

    df_input = pd.read_csv(csv_path)
    records = []

    for _, row in df_input.iterrows():
        name = str(row["name"])
        smiles = str(row["smiles"])
        logger.info("Processing ligand: %s", name)

        # --- Parse SMILES ---
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Invalid SMILES for %s — skipping: %s", name, smiles)
            continue

        # Add explicit hydrogens before 3D embedding
        # Hydrogens are needed for accurate docking (H-bond donors, charges)
        mol = AllChem.AddHs(mol)

        # --- 3D conformer generation and minimisation ---
        mol = _embed_and_minimise(mol)
        if mol is None:
            logger.warning("3D embedding failed for %s — skipping", name)
            continue

        # --- Write SDF (canonical 3D format for structure exchange) ---
        sdf_path = ligands_dir / f"{name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(mol)
        writer.close()

        # --- Write PDB (human-readable; useful for visual inspection) ---
        pdb_path = ligands_dir / f"{name}.pdb"
        Chem.MolToPDBFile(mol, str(pdb_path))

        # --- Generate PDBQT via Meeko ---
        pdbqt_path = ligands_dir / f"{name}.pdbqt"
        meeko_ok = _run_meeko(sdf_path, pdbqt_path)
        if not meeko_ok:
            logger.warning(
                "PDBQT generation failed for %s — ligand will not be dockable", name
            )

        # --- Compute descriptors ---
        # Remove Hs for descriptor calculation (standard practice)
        mol_noH = AllChem.RemoveHs(mol)
        desc = _compute_descriptors(mol_noH, name)
        desc["pdbqt_ready"] = meeko_ok
        records.append(desc)

        logger.info(
            "%s: MW=%.1f, logP=%.2f, Lipinski=%s, Veber=%s",
            name, desc["MW"], desc["MolLogP"],
            desc["Lipinski_pass"], desc["Veber_pass"],
        )

    if not records:
        logger.error("No ligands successfully processed")
        return pd.DataFrame()

    df_out = pd.DataFrame(records)
    out_path = ligands_dir / "ligand_properties.csv"
    df_out.to_csv(out_path, index=False)
    logger.info("Ligand properties written to %s", out_path)

    return df_out
