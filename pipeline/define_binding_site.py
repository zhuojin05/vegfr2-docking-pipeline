"""
Stage 3: Extract the VEGFR2 binding site from the sorafenib co-crystal structure.

Uses the 3WZE co-crystal complex (VEGFR2 + sorafenib) to automatically define
the docking box. Sorafenib's HETATM residue name in 3WZE is 'BAX' (residue 1201,
chain A) — verified by grepping the raw PDB. Do NOT assume SFN/SOR/SF4.

The box is centred on the ligand centroid and padded by config['docking']['padding']
on each side. This approach is standard for structure-based docking when a
co-crystal ligand is available.

Note: the ICL practical manual quotes a centroid of (~1.3, ~6.8, ~9.1 Å), but
that value is in the 4ASD coordinate frame. The pipeline extracts the box from
3WZE, which has its own frame — the values are not comparable.
"""

import json
import logging
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer

logger = logging.getLogger(__name__)


def _find_ligand_residue(structure, residue_name: str) -> list[np.ndarray]:
    """
    Search all models/chains for a HETATM residue with the given name.

    BioPython structure hierarchy: Structure → Model → Chain → Residue → Atom.
    HETATM residues are stored with residue id flag 'H_<name>' in BioPython.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
    residue_name : str
        Three-letter HETATM residue name to search for (e.g. 'BAX').

    Returns
    -------
    list of np.ndarray
        Heavy atom coordinates (excludes hydrogen atoms).
    """
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname.strip() == residue_name:
                    logger.info(
                        "Found residue %s in chain %s, residue id %s",
                        residue_name, chain.id, residue.id
                    )
                    for atom in residue:
                        # Skip hydrogen atoms — only heavy atoms define the
                        # pharmacophoric volume of the ligand
                        if atom.element != "H" and not atom.name.startswith("H"):
                            coords.append(atom.coord)
    return coords


def _get_superimposer(receptor_raw_pdb: Path, reference_pdb: Path) -> Superimposer:
    """
    Compute the rigid-body transformation that maps 3WZE coordinates into the
    4ASD coordinate frame, using shared Cα atoms (Kabsch algorithm).

    Uses the raw receptor PDB (original residue numbering 807–1168), not the
    PDBFixer-cleaned version (renumbered 1–305), so residue numbers overlap
    with 3WZE chain A residues 815–1167. Atom coordinates are identical.

    BioPython convention: sup.rotran = (R, t), row-vector form.
    Transform: coord_4asd = coord_3wze @ R + t
    """
    parser = PDBParser(QUIET=True)
    receptor = parser.get_structure("receptor", str(receptor_raw_pdb))
    reference = parser.get_structure("reference", str(reference_pdb))

    def _ca_map(structure):
        result = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        result[(chain.id, residue.id[1])] = residue["CA"]
        return result

    ca_rec = _ca_map(receptor)    # 4ASD: chain A, residues 807–1168
    ca_ref = _ca_map(reference)   # 3WZE: chain A, residues 815–1167
    common = sorted(set(ca_rec) & set(ca_ref))

    if len(common) < 30:
        raise ValueError(
            f"Only {len(common)} common Cα residues between "
            f"{receptor_raw_pdb.name} and {reference_pdb.name}. "
            "Check that both PDB files are the same protein."
        )

    sup = Superimposer()
    sup.set_atoms([ca_rec[k] for k in common], [ca_ref[k] for k in common])

    logger.info(
        "Superimposed %s onto %s: %d Cα pairs, RMSD = %.3f Å",
        reference_pdb.name, receptor_raw_pdb.name, len(common), sup.rms
    )
    if sup.rms > 3.0:
        logger.warning(
            "Superimposition RMSD %.3f Å is high — binding site may be inaccurate.",
            sup.rms
        )
    return sup


def define_binding_site(config: dict) -> dict:
    """
    Stage 3: Compute docking box from the sorafenib (BAX) co-crystal pose in 3WZE.

    Parameters
    ----------
    config : dict
        Pipeline config loaded from config.yaml.

    Returns
    -------
    dict
        Binding site parameters: center_x/y/z, size_x/y/z, padding, ligand_residue,
        reference_pdb. Also written to data/prepared/binding_site.json and
        data/prepared/config_template.txt.
    """
    raw_dir = Path(config["paths"]["raw_data"])
    prepared_dir = Path(config["paths"]["prepared"])
    prepared_dir.mkdir(parents=True, exist_ok=True)

    reference_pdb = config["receptor"]["reference_complex"]
    pdb_path = raw_dir / reference_pdb
    if not pdb_path.exists():
        raise FileNotFoundError(f"Reference complex PDB not found: {pdb_path}")

    padding = float(config["docking"]["padding"])

    # Sorafenib's actual HETATM name in 3WZE — verified by:
    #   grep HETATM data/raw/3WZE.pdb | grep -v HOH
    # Result: residue BAX, chain A, residue number 1201
    ligand_name = "BAX"

    logger.info("Loading reference complex: %s", pdb_path)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("3WZE", str(pdb_path))

    coords = _find_ligand_residue(structure, ligand_name)

    if not coords:
        # Fallback: try common sorafenib names in case PDB was re-deposited
        fallback_names = ["SFN", "SOR", "SF4", "SRF"]
        for alt_name in fallback_names:
            coords = _find_ligand_residue(structure, alt_name)
            if coords:
                ligand_name = alt_name
                logger.warning(
                    "BAX not found; using fallback residue name '%s'", alt_name
                )
                break

    if not coords:
        raise ValueError(
            f"Could not find sorafenib ligand (tried BAX, SFN, SOR, SF4, SRF) "
            f"in {pdb_path}. Run: grep HETATM data/raw/3WZE.pdb to inspect."
        )

    coords_arr_3wze = np.array(coords)
    logger.info(
        "Found %d heavy atoms for %s in 3WZE frame (centroid: %.3f, %.3f, %.3f)",
        len(coords_arr_3wze), ligand_name, *coords_arr_3wze.mean(axis=0)
    )

    # Transform from 3WZE coordinate frame into 4ASD frame.
    # These are different crystal structures with different unit cells —
    # direct coordinate comparison is meaningless without alignment.
    receptor_raw_pdb = raw_dir / config["receptor"]["pdb"]
    sup = _get_superimposer(receptor_raw_pdb, pdb_path)
    R, t = sup.rotran
    coords_arr = coords_arr_3wze @ R + t   # now in 4ASD frame

    logger.info(
        "BAX centroid in 4ASD frame: x=%.3f, y=%.3f, z=%.3f", *coords_arr.mean(axis=0)
    )

    # Centroid = geometric centre of the co-crystal ligand (in 4ASD frame)
    centroid = coords_arr.mean(axis=0)

    # Box size = ligand extent + padding on each side
    # The padding ensures the full binding pocket is sampled, not just the
    # footprint of the co-crystal pose
    min_coords = coords_arr.min(axis=0)
    max_coords = coords_arr.max(axis=0)
    box_size = (max_coords - min_coords) + 2 * padding

    site = {
        "center_x": float(round(centroid[0], 3)),
        "center_y": float(round(centroid[1], 3)),
        "center_z": float(round(centroid[2], 3)),
        "size_x": float(round(box_size[0], 3)),
        "size_y": float(round(box_size[1], 3)),
        "size_z": float(round(box_size[2], 3)),
        "padding": padding,
        "ligand_residue": ligand_name,
        "reference_pdb": reference_pdb,
        "coordinate_frame": config["receptor"]["pdb"],
        "superimposition_rmsd_angstrom": float(round(sup.rms, 4)),
    }

    # Write JSON (machine-readable, used by run_docking.py)
    json_path = prepared_dir / "binding_site.json"
    with open(json_path, "w") as f:
        json.dump(site, f, indent=2)
    logger.info("Binding site parameters written to %s", json_path)

    # Write VINA-format config template (human-readable reference)
    txt_path = prepared_dir / "config_template.txt"
    with open(txt_path, "w") as f:
        f.write(f"center_x = {site['center_x']}\n")
        f.write(f"center_y = {site['center_y']}\n")
        f.write(f"center_z = {site['center_z']}\n")
        f.write(f"size_x = {site['size_x']}\n")
        f.write(f"size_y = {site['size_y']}\n")
        f.write(f"size_z = {site['size_z']}\n")
    logger.info("VINA config template written to %s", txt_path)

    logger.info(
        "Binding site summary:\n"
        "  Centre: (%.3f, %.3f, %.3f)\n"
        "  Box:    (%.1f × %.1f × %.1f) Å",
        site["center_x"], site["center_y"], site["center_z"],
        site["size_x"], site["size_y"], site["size_z"],
    )

    return site
