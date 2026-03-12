"""
Tests for Stage 2: prepare_receptor.py

Tests PDBFixer loading and heteroatom removal on 4ASD.pdb.
Does not test Meeko PDBQT generation (requires external binary).
"""

import pytest
from pathlib import Path

# Resolve project root relative to this test file
PROJECT_ROOT = Path(__file__).parent.parent
RAW_PDB = PROJECT_ROOT / "data" / "raw" / "4ASD.pdb"


def _pdbfixer_available() -> bool:
    try:
        import pdbfixer  # noqa: F401
        import openmm  # noqa: F401
        return True
    except ImportError:
        return False


PDBFIXER_AVAILABLE = _pdbfixer_available()
SKIP_MSG = "pdbfixer/openmm not installed (conda-forge required on arm64)"


@pytest.mark.skipif(not RAW_PDB.exists(), reason="4ASD.pdb not found in data/raw/")
@pytest.mark.skipif(not PDBFIXER_AVAILABLE, reason=SKIP_MSG)
class TestPDBFixerLoads:
    def test_loads_raw_pdb(self):
        """PDBFixer should load 4ASD.pdb without error."""
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(RAW_PDB))
        assert fixer.topology is not None

    def test_chain_count(self):
        """4ASD has at least one protein chain."""
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(RAW_PDB))
        chains = list(fixer.topology.chains())
        assert len(chains) >= 1

    def test_atom_count_nonzero(self):
        """Structure should have a non-trivial number of atoms."""
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(RAW_PDB))
        atoms = list(fixer.topology.atoms())
        # 4ASD is ~320 residues; expect thousands of atoms
        assert len(atoms) > 1000

    def test_missing_residues_runs(self):
        """findMissingResidues should run without error."""
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(RAW_PDB))
        fixer.findMissingResidues()  # should not raise

    def test_find_and_add_missing_atoms(self):
        """findMissingAtoms and addMissingAtoms should not crash."""
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(RAW_PDB))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()


@pytest.mark.skipif(not RAW_PDB.exists(), reason="4ASD.pdb not found in data/raw/")
@pytest.mark.skipif(not PDBFIXER_AVAILABLE, reason=SKIP_MSG)
class TestRemoveHeterogens:
    def test_remove_heterogens(self):
        """After removeHeterogens, topology should have no heteroatom residues."""
        from pdbfixer import PDBFixer
        import openmm.app as app
        import io

        fixer = PDBFixer(filename=str(RAW_PDB))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.4)
        fixer.removeHeterogens(keepWater=False)

        # Write to a string buffer and check for HETATM records
        buf = io.StringIO()
        app.PDBFile.writeFile(fixer.topology, fixer.positions, buf)
        pdb_text = buf.getvalue()

        # No HETATM records should remain (waters and ligands removed)
        hetatm_lines = [l for l in pdb_text.splitlines() if l.startswith("HETATM")]
        assert len(hetatm_lines) == 0, (
            f"Expected 0 HETATM records after removal, found {len(hetatm_lines)}"
        )

    def test_hydrogens_added(self):
        """After addMissingHydrogens, atom count should increase."""
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=str(RAW_PDB))
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        n_before = len(list(fixer.topology.atoms()))
        fixer.addMissingHydrogens(pH=7.4)
        n_after = len(list(fixer.topology.atoms()))
        assert n_after > n_before, "addMissingHydrogens should increase atom count"
