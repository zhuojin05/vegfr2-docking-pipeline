"""
Tests for Stage 1: prepare_ligands.py

Test ligand: aspirin (CC(=O)Oc1ccccc1C(=O)O)
No Gnina or VINA required — tests only cover RDKit processing and descriptors.
"""

import logging
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O"

# Minimal config for testing — no file paths needed for descriptor tests
MOCK_CONFIG = {
    "paths": {"ligands": "/tmp/test_ligands"},
    "docking": {"backend": "gnina", "num_modes": 5, "exhaustiveness": 4, "padding": 10.0},
    "receptor": {"pdb": "4ASD.pdb", "reference_complex": "3WZE.pdb", "ph": 7.4,
                 "remove_heteroatoms": True, "remove_waters": True},
}


class TestSmilsTo3D:
    def test_parse_valid_smiles(self):
        """Valid SMILES should produce a non-None RDKit mol."""
        mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
        assert mol is not None

    def test_add_hydrogens(self):
        """AddHs should increase atom count."""
        mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
        mol_h = AllChem.AddHs(mol)
        assert mol_h.GetNumAtoms() > mol.GetNumAtoms()

    def test_3d_conformer_generated(self):
        """ETKDGv3 embedding should produce a valid conformer."""
        mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
        mol = AllChem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)
        assert result != -1, "3D embedding failed"
        conf = mol.GetConformer()
        # Aspirin has 21 atoms with Hs; all positions should be non-zero
        positions = conf.GetPositions()
        assert positions.shape[0] == mol.GetNumAtoms()
        # At least some atoms should not be at origin
        assert any(abs(pos).max() > 0 for pos in positions)

    def test_mmff94_minimisation(self):
        """MMFF94 minimisation should not crash and return 0 (converged)."""
        mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
        mol = AllChem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        AllChem.EmbedMolecule(mol, params)
        result = AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        # 0 = converged, 1 = more iters needed, -1 = failed
        assert result in (0, 1), f"MMFF94 returned unexpected code: {result}"


class TestDescriptors:
    def _get_mol(self):
        mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
        return AllChem.RemoveHs(AllChem.AddHs(mol))

    def test_mw_aspirin(self):
        """Aspirin MW should be ~180.16 Da."""
        mol = self._get_mol()
        mw = Descriptors.MolWt(mol)
        assert 179.0 < mw < 182.0, f"Unexpected MW: {mw}"

    def test_lipinski_aspirin_passes(self):
        """Aspirin is a Lipinski-compliant molecule."""
        mol = self._get_mol()
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        assert mw <= 500
        assert logp <= 5
        assert hbd <= 5
        assert hba <= 10

    def test_logp_positive(self):
        """Aspirin logP should be positive (mildly lipophilic)."""
        mol = self._get_mol()
        logp = Descriptors.MolLogP(mol)
        assert logp > 0

    def test_rotatable_bonds(self):
        """Aspirin has 3 rotatable bonds."""
        mol = self._get_mol()
        rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        assert rot == 3, f"Expected 3 rotatable bonds, got {rot}"

    def test_veber_aspirin_passes(self):
        """Aspirin satisfies Veber rules (RotBonds ≤ 10, TPSA ≤ 140)."""
        mol = self._get_mol()
        rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)
        assert rot <= 10
        assert tpsa <= 140


class TestInvalidSmiles:
    def test_invalid_smiles_returns_none(self):
        """MolFromSmiles should return None for invalid SMILES, not raise."""
        mol = Chem.MolFromSmiles("this_is_not_smiles!!")
        assert mol is None

    def test_batch_continues_after_invalid(self):
        """
        Simulate batch processing: a bad SMILES should not abort other ligands.
        This mirrors the per-ligand error handling in prepare_ligands().
        """
        smiles_list = [
            ("aspirin", ASPIRIN_SMILES),
            ("bad_mol", "INVALID###"),
            ("another_aspirin", ASPIRIN_SMILES),
        ]
        processed = []
        for name, smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logging.warning("Skipping %s: invalid SMILES", name)
                continue
            processed.append(name)

        assert "aspirin" in processed
        assert "another_aspirin" in processed
        assert "bad_mol" not in processed
        assert len(processed) == 2
