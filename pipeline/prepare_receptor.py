"""
Stage 2: Prepare the VEGFR2 receptor (4ASD) for docking.

Uses PDBFixer to:
  - Fill missing residues and atoms
  - Add hydrogens at physiological pH 7.4
  - Remove HETATM records (co-crystal ligands, buffer molecules) and waters

Then calls Meeko's mk_prepare_receptor.py to produce the PDBQT file required
by AutoDock VINA and Gnina.

pdbfixer and openmm are installed via uv sync (PyPI arm64 wheels available
from the OpenMM project since v8.x).
"""

import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def prepare_receptor(config: dict) -> Path:
    """
    Stage 2: Fix and protonate the VEGFR2 receptor; generate PDBQT.

    Parameters
    ----------
    config : dict
        Pipeline config loaded from config.yaml.

    Returns
    -------
    Path
        Path to the cleaned PDB file (data/prepared/4ASD_clean.pdb).
        The PDBQT is written alongside it.
    """
    try:
        from pdbfixer import PDBFixer
        import openmm.app as app
    except ImportError as exc:
        raise ImportError(
            "pdbfixer/openmm not found. Install with: uv sync"
        ) from exc

    raw_dir = Path(config["paths"]["raw_data"])
    prepared_dir = Path(config["paths"]["prepared"])
    prepared_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = config["receptor"]["pdb"]
    ph = float(config["receptor"]["ph"])

    raw_path = raw_dir / receptor_pdb
    if not raw_path.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {raw_path}")

    logger.info("Loading receptor: %s", raw_path)
    fixer = PDBFixer(filename=str(raw_path))

    # Identify and fill structural gaps
    # 4ASD is a high-resolution structure (2.0 Å) but may have disordered loops
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    logger.info("Missing residues/atoms identified and added")

    # Add hydrogens at physiological pH 7.4
    # Protonation state affects H-bond donor/acceptor assignment in docking
    # Histidine protonation is particularly important near the ATP-binding hinge
    fixer.addMissingHydrogens(pH=ph)
    logger.info("Hydrogens added at pH %.1f", ph)

    # Remove all HETATM records: co-crystal ligands, ions, buffer molecules
    # keepWater=False also removes crystallographic waters
    # We want the apo receptor for unbiased docking
    fixer.removeHeterogens(keepWater=False)
    logger.info("Heterogens and waters removed")

    # Write cleaned PDB using OpenMM's writer
    stem = receptor_pdb.replace(".pdb", "")
    clean_pdb_path = prepared_dir / f"{stem}_clean.pdb"
    with open(clean_pdb_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    logger.info("Cleaned receptor written to %s", clean_pdb_path)

    # Generate PDBQT via Meeko
    # mk_prepare_receptor.py handles AutoDock atom typing and charge assignment
    # In Meeko ≥0.5, -o sets the output *basename* for JSON/GPF files.
    # The -p flag is required to trigger PDBQT output; without it, no .pdbqt
    # is written even when the command exits 0.
    pdbqt_path = prepared_dir / f"{stem}.pdbqt"
    pdbqt_stem = prepared_dir / stem   # Meeko appends .pdbqt itself when -p is used
    cmd = [
        "mk_prepare_receptor.py",
        "-i", str(clean_pdb_path),
        "-o", str(pdbqt_stem),         # basename; Meeko appends .pdbqt
        "-p",                          # write PDBQT output (required in Meeko ≥0.5)
    ]
    logger.info("Running Meeko receptor preparation: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(
            "mk_prepare_receptor.py failed:\n%s", result.stderr.strip()
        )
        raise RuntimeError(
            f"Receptor PDBQT generation failed. stderr:\n{result.stderr}"
        )

    # Verify the file was actually created; Meeko exits 0 even if nothing is
    # written in some edge cases (e.g., unsupported input), so check explicitly.
    if not pdbqt_path.exists():
        raise RuntimeError(
            f"mk_prepare_receptor.py completed without error but did not write "
            f"{pdbqt_path}. Check that the -p flag is supported by your Meeko version."
        )
    logger.info("Receptor PDBQT written to %s", pdbqt_path)
    return clean_pdb_path
