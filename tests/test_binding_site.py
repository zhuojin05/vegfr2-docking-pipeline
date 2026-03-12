"""
Tests for Stage 3: define_binding_site.py

Verifies that the BAX sorafenib residue can be found in 3WZE.pdb and that
the computed centroid and JSON output are correct.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser

PROJECT_ROOT = Path(__file__).parent.parent
WZE_PDB = PROJECT_ROOT / "data" / "raw" / "3WZE.pdb"



@pytest.mark.skipif(not WZE_PDB.exists(), reason="3WZE.pdb not found in data/raw/")
class TestBAXResidue:
    def _load_structure(self):
        parser = PDBParser(QUIET=True)
        return parser.get_structure("3WZE", str(WZE_PDB))

    def test_bax_residue_found(self):
        """3WZE.pdb must contain a HETATM residue named BAX (sorafenib)."""
        structure = self._load_structure()
        found = False
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname.strip() == "BAX":
                        found = True
                        break
        assert found, (
            "Residue 'BAX' not found in 3WZE.pdb. "
            "Run: grep HETATM data/raw/3WZE.pdb to inspect HETATM names."
        )

    def test_bax_heavy_atoms(self):
        """BAX residue should have multiple heavy atoms (sorafenib has 27)."""
        structure = self._load_structure()
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname.strip() == "BAX":
                        heavy_atoms = [
                            a for a in residue
                            if a.element != "H" and not a.name.startswith("H")
                        ]
                        assert len(heavy_atoms) > 10, (
                            f"Expected >10 heavy atoms in BAX, found {len(heavy_atoms)}"
                        )
                        return
        pytest.fail("BAX not found")

    def test_centroid_in_range(self):
        """BAX centroid should be within the 3WZE structure (sanity: non-zero, finite)."""
        structure = self._load_structure()
        coords = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.resname.strip() == "BAX":
                        for atom in residue:
                            if atom.element != "H" and not atom.name.startswith("H"):
                                coords.append(atom.coord)
        assert len(coords) > 0
        centroid = np.mean(coords, axis=0)
        assert np.all(np.isfinite(centroid)), "BAX centroid contains non-finite values"
        assert np.linalg.norm(centroid) > 1.0, "BAX centroid suspiciously close to origin"


@pytest.mark.skipif(not WZE_PDB.exists(), reason="3WZE.pdb not found in data/raw/")
class TestJSONOutput:
    def _mock_config(self, prepared_dir: str) -> dict:
        return {
            "paths": {
                "raw_data": str(PROJECT_ROOT / "data" / "raw"),
                "prepared": prepared_dir,
                "ligands": str(PROJECT_ROOT / "data" / "ligands"),
                "results": str(PROJECT_ROOT / "data" / "results"),
                "logs": str(PROJECT_ROOT / "logs"),
            },
            "docking": {
                "backend": "gnina",
                "num_modes": 10,
                "exhaustiveness": 8,
                "padding": 10.0,
            },
            "receptor": {
                "pdb": "4ASD.pdb",
                "reference_complex": "3WZE.pdb",
                "ph": 7.4,
                "remove_heteroatoms": True,
                "remove_waters": True,
            },
        }

    def test_json_output_keys(self):
        """define_binding_site() should write a JSON with the expected keys."""
        from pipeline.define_binding_site import define_binding_site

        with tempfile.TemporaryDirectory() as tmpdir:
            config = self._mock_config(tmpdir)
            site = define_binding_site(config)

            required_keys = [
                "center_x", "center_y", "center_z",
                "size_x", "size_y", "size_z",
                "padding", "ligand_residue", "reference_pdb",
                "coordinate_frame", "superimposition_rmsd_angstrom",
            ]
            for key in required_keys:
                assert key in site, f"Missing key in binding site dict: {key}"

            # Verify JSON file was written
            json_path = Path(tmpdir) / "binding_site.json"
            assert json_path.exists(), "binding_site.json not written"

            with open(json_path) as f:
                loaded = json.load(f)
            for key in required_keys:
                assert key in loaded, f"Missing key in binding_site.json: {key}"

    def test_json_values_are_finite(self):
        """All numeric binding site parameters should be finite."""
        from pipeline.define_binding_site import define_binding_site

        with tempfile.TemporaryDirectory() as tmpdir:
            config = self._mock_config(tmpdir)
            site = define_binding_site(config)

            for key in ["center_x", "center_y", "center_z",
                        "size_x", "size_y", "size_z"]:
                assert isinstance(site[key], float), f"{key} should be float"
                assert np.isfinite(site[key]), f"{key} is not finite: {site[key]}"

    def test_box_size_positive(self):
        """Box dimensions must be positive (padding > 0 ensures this)."""
        from pipeline.define_binding_site import define_binding_site

        with tempfile.TemporaryDirectory() as tmpdir:
            config = self._mock_config(tmpdir)
            site = define_binding_site(config)

            assert site["size_x"] > 0
            assert site["size_y"] > 0
            assert site["size_z"] > 0

    def test_config_template_written(self):
        """config_template.txt should be written alongside binding_site.json."""
        from pipeline.define_binding_site import define_binding_site

        with tempfile.TemporaryDirectory() as tmpdir:
            config = self._mock_config(tmpdir)
            define_binding_site(config)

            txt_path = Path(tmpdir) / "config_template.txt"
            assert txt_path.exists(), "config_template.txt not written"

            content = txt_path.read_text()
            for field in ["center_x", "center_y", "center_z",
                          "size_x", "size_y", "size_z"]:
                assert field in content, f"Field '{field}' missing from config_template.txt"
