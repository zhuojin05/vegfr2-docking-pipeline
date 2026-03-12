"""
Stage 4: Run molecular docking with Gnina (preferred) or AutoDock VINA.

Gnina uses a CNN-based scoring function trained on PDBbind, which is more
accurate than VINA's empirical scoring for kinase targets like VEGFR2.
VINA is provided as a comparison backend (faster on CPU, empirical scoring).

Output scores:
  - Gnina: CNNscore (CNN confidence 0–1), CNNaffinity (predicted pKd),
    minimizedAffinity (kcal/mol, same scale as VINA)
  - VINA: affinity (kcal/mol), more negative = better binding

APPLE SILICON NOTE:
  Gnina requires CUDA 12.0+, which is unavailable on macOS. There is no
  native arm64 binary and no Homebrew formula. The only working option on
  Apple Silicon is the 'gnina_docker' backend:

    docker pull gnina/gnina:latest
    # then set: docking.backend: gnina_docker  in config.yaml

  Docker Desktop emulates x86_64 Linux via Rosetta 2. Inside the container,
  gnina's --no_gpu flag disables CUDA initialisation entirely, so it runs as
  a CPU-only process. Expect ~10–30 min per ligand (x86_64 emulation + CNN
  scoring). For faster local runs, use 'vina' (default for Apple Silicon).
"""

import json
import logging
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Gnina backend
# ---------------------------------------------------------------------------

def _parse_gnina_pdbqt(pdbqt_path: Path) -> pd.DataFrame:
    """
    Parse Gnina output PDBQT to extract per-pose scores.

    Gnina annotates each MODEL block with REMARK lines:
      REMARK VINA RESULT:   <minimizedAffinity>  <rmsd_lb>  <rmsd_ub>
      REMARK CNNscore <value>
      REMARK CNNaffinity <value>

    Parameters
    ----------
    pdbqt_path : Path

    Returns
    -------
    pd.DataFrame
        Columns: pose, minimized_affinity, rmsd_lb, rmsd_ub, cnn_score, cnn_affinity.
    """
    records = []
    current_pose = 0
    current = {}

    with open(pdbqt_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("MODEL"):
                current_pose += 1
                current = {"pose": current_pose}
            elif line.startswith("REMARK VINA RESULT:"):
                parts = line.split()[3:]
                if len(parts) >= 3:
                    current["minimized_affinity"] = float(parts[0])
                    current["rmsd_lb"] = float(parts[1])
                    current["rmsd_ub"] = float(parts[2])
            elif line.startswith("REMARK CNNscore"):
                m = re.search(r"CNNscore\s+([\d.eE+\-]+)", line)
                if m:
                    current["cnn_score"] = float(m.group(1))
            elif line.startswith("REMARK CNNaffinity"):
                m = re.search(r"CNNaffinity\s+([\d.eE+\-]+)", line)
                if m:
                    current["cnn_affinity"] = float(m.group(1))
            elif line.startswith("ENDMDL") and current:
                records.append(current)
                current = {}

    return pd.DataFrame(records)


def _run_gnina(ligand_name: str, site: dict, config: dict) -> pd.DataFrame:
    """
    Run Gnina docking for one ligand and return scored poses as a DataFrame.

    Parameters
    ----------
    ligand_name : str
    site : dict
        Binding site parameters from define_binding_site().
    config : dict

    Returns
    -------
    pd.DataFrame
    """
    if shutil.which("gnina") is None:
        raise RuntimeError(
            "gnina not found in PATH. Install options:\n"
            "  brew install gnina          # if arm64 binary available\n"
            "  Download from https://github.com/gnina/gnina/releases\n"
            "  arch -x86_64 brew install gnina  # Rosetta 2 fallback\n"
            "See README.md for details."
        )

    prepared_dir = Path(config["paths"]["prepared"])
    ligands_dir = Path(config["paths"]["ligands"])
    results_dir = Path(config["paths"]["results"])
    results_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = prepared_dir / "4ASD.pdbqt"
    ligand_pdbqt = ligands_dir / f"{ligand_name}.pdbqt"

    for p in [receptor_pdbqt, ligand_pdbqt]:
        if not p.exists():
            raise FileNotFoundError(f"Required file not found: {p}")

    out_pdbqt = results_dir / f"{ligand_name}_gnina_poses.pdbqt"
    out_log = results_dir / f"{ligand_name}_gnina.log"

    cmd = [
        "gnina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_pdbqt),
        "--center_x", str(site["center_x"]),
        "--center_y", str(site["center_y"]),
        "--center_z", str(site["center_z"]),
        "--size_x", str(site["size_x"]),
        "--size_y", str(site["size_y"]),
        "--size_z", str(site["size_z"]),
        "--num_modes", str(config["docking"]["num_modes"]),
        "--exhaustiveness", str(config["docking"]["exhaustiveness"]),
        "--scoring", "cnnscore",
        "--out", str(out_pdbqt),
        "--log", str(out_log),
    ]

    logger.info("Running Gnina for %s (CPU mode — may be slow)", ligand_name)
    logger.debug("Command: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("Gnina failed for %s:\n%s", ligand_name, result.stderr[:500])
        raise RuntimeError(f"Gnina docking failed for {ligand_name}")

    if not out_pdbqt.exists():
        raise FileNotFoundError(f"Gnina output not created: {out_pdbqt}")

    df = _parse_gnina_pdbqt(out_pdbqt)
    df["ligand"] = ligand_name
    df["backend"] = "gnina"
    # Rename for unified schema
    df["score"] = df.get("cnn_score", float("nan"))

    return df


# ---------------------------------------------------------------------------
# Gnina-via-Docker backend (Apple Silicon workaround)
# ---------------------------------------------------------------------------

def _run_gnina_docker(ligand_name: str, site: dict, config: dict) -> pd.DataFrame:
    """
    Run Gnina docking inside a Docker container (x86_64 Rosetta 2 emulation).

    This is the supported path for Apple Silicon, where Gnina cannot run natively
    due to the CUDA 12.0+ requirement. Docker Desktop emulates x86_64 Linux via
    Rosetta 2; gnina's --no_gpu flag disables CUDA initialisation entirely.

    Prerequisites (one-time manual setup):
      1. Install Docker Desktop for Mac (arm64):
         https://www.docker.com/products/docker-desktop/
      2. Pull the gnina image:
         docker pull gnina/gnina:latest

    Parameters
    ----------
    ligand_name : str
    site : dict
        Binding site parameters from define_binding_site().
    config : dict

    Returns
    -------
    pd.DataFrame
        Same schema as _run_gnina() (CNNscore, CNNaffinity, minimizedAffinity).
    """
    if shutil.which("docker") is None:
        raise RuntimeError(
            "Docker not found in PATH. Install Docker Desktop for Mac (arm64):\n"
            "  https://www.docker.com/products/docker-desktop/\n"
            "Then pull the gnina image:\n"
            "  docker pull gnina/gnina:latest\n"
            "Finally set 'docking.backend: gnina_docker' in config.yaml."
        )

    # Resolve project root (two levels up from this file: pipeline/ → project root)
    project_root = Path(__file__).resolve().parents[1]

    prepared_dir = Path(config["paths"]["prepared"])
    ligands_dir = Path(config["paths"]["ligands"])
    results_dir = Path(config["paths"]["results"])
    results_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = prepared_dir / "4ASD.pdbqt"
    ligand_pdbqt = ligands_dir / f"{ligand_name}.pdbqt"

    for p in [receptor_pdbqt, ligand_pdbqt]:
        if not p.exists():
            raise FileNotFoundError(f"Required file not found: {p}")

    out_pdbqt = results_dir / f"{ligand_name}_gnina_poses.pdbqt"
    out_log = results_dir / f"{ligand_name}_gnina.log"

    # Paths inside the container are relative to /workspace (= project_root on host)
    receptor_in_container = "/workspace" / receptor_pdbqt.resolve().relative_to(project_root)
    ligand_in_container = "/workspace" / ligand_pdbqt.resolve().relative_to(project_root)
    out_in_container = "/workspace" / out_pdbqt.resolve().relative_to(project_root)
    log_in_container = "/workspace" / out_log.resolve().relative_to(project_root)

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{project_root}:/workspace",
        "gnina/gnina:latest",
        "gnina",
        "--no_gpu",  # disable CUDA initialisation — required in Docker on macOS
        "--receptor", str(receptor_in_container),
        "--ligand", str(ligand_in_container),
        "--center_x", str(site["center_x"]),
        "--center_y", str(site["center_y"]),
        "--center_z", str(site["center_z"]),
        "--size_x", str(site["size_x"]),
        "--size_y", str(site["size_y"]),
        "--size_z", str(site["size_z"]),
        "--num_modes", str(config["docking"]["num_modes"]),
        "--exhaustiveness", str(config["docking"]["exhaustiveness"]),
        "--scoring", "cnnscore",
        "--out", str(out_in_container),
        "--log", str(log_in_container),
    ]

    # x86_64 emulation via Rosetta 2 + CPU CNN scoring is inherently slow
    logger.warning(
        "Running gnina via Docker (x86_64 Rosetta 2 emulation — expect ~10–30 min per ligand)"
    )
    logger.info("Running Gnina (Docker) for %s", ligand_name)
    logger.debug("Command: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("Gnina (Docker) failed for %s:\n%s", ligand_name, result.stderr[:500])
        raise RuntimeError(f"Gnina (Docker) docking failed for {ligand_name}")

    if not out_pdbqt.exists():
        raise FileNotFoundError(f"Gnina Docker output not created: {out_pdbqt}")

    df = _parse_gnina_pdbqt(out_pdbqt)
    df["ligand"] = ligand_name
    df["backend"] = "gnina_docker"
    df["score"] = df.get("cnn_score", float("nan"))

    return df


# ---------------------------------------------------------------------------
# VINA backend
# ---------------------------------------------------------------------------

def _parse_vina_rmsd(pdbqt_path: Path) -> list[tuple[float, float]]:
    """
    Extract rmsd_lb and rmsd_ub from a multi-model VINA PDBQT file.

    VINA writes one MODEL/ENDMDL block per pose; each block contains:
        REMARK VINA RESULT:  affinity  rmsd_lb  rmsd_ub
    The Python API energies() does not expose RMSD, so we read it here.
    """
    results = []
    with open(pdbqt_path) as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                # parts: ['REMARK', 'VINA', 'RESULT:', affinity, rmsd_lb, rmsd_ub]
                results.append((float(parts[4]), float(parts[5])))
    return results


def _run_vina(ligand_name: str, site: dict, config: dict) -> pd.DataFrame:
    """
    Run AutoDock VINA docking via Python bindings.

    VINA uses an empirical scoring function (kcal/mol); more negative = tighter
    predicted binding. It is faster than Gnina on CPU and useful as a comparison.

    Parameters
    ----------
    ligand_name : str
    site : dict
    config : dict

    Returns
    -------
    pd.DataFrame
    """
    try:
        from vina import Vina
    except ImportError as exc:
        raise ImportError(
            "AutoDock VINA Python bindings not found.\n"
            "Install: pip install vina --no-binary vina"
        ) from exc

    prepared_dir = Path(config["paths"]["prepared"])
    ligands_dir = Path(config["paths"]["ligands"])
    results_dir = Path(config["paths"]["results"])
    results_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdbqt = prepared_dir / "4ASD.pdbqt"
    ligand_pdbqt = ligands_dir / f"{ligand_name}.pdbqt"

    for p in [receptor_pdbqt, ligand_pdbqt]:
        if not p.exists():
            raise FileNotFoundError(f"Required file not found: {p}")

    v = Vina(sf_name="vina")
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))

    center = [site["center_x"], site["center_y"], site["center_z"]]
    box_size = [site["size_x"], site["size_y"], site["size_z"]]
    v.compute_vina_maps(center=center, box_size=box_size)

    n_poses = int(config["docking"]["num_modes"])
    exhaustiveness = int(config["docking"]["exhaustiveness"])

    logger.info("Running VINA for %s", ligand_name)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    out_pdbqt = results_dir / f"{ligand_name}_vina_poses.pdbqt"
    v.write_poses(str(out_pdbqt), n_poses=n_poses, overwrite=True)

    # energies() returns array of shape (n_poses, 6):
    # [total, inter, intra, torsion, intra_best_pose, total_intra]
    energies = v.energies(n_poses=n_poses)

    # RMSD is not returned by energies(); parse it from the PDBQT we just wrote.
    # Pose 1 always has rmsd_lb = rmsd_ub = 0.0 (reference pose by definition).
    rmsd_vals = _parse_vina_rmsd(out_pdbqt)

    records = []
    for i, row in enumerate(energies):
        rmsd_lb, rmsd_ub = rmsd_vals[i] if i < len(rmsd_vals) else (float("nan"), float("nan"))
        records.append({
            "pose": i + 1,
            "ligand": ligand_name,
            "backend": "vina",
            "score": float(row[0]),         # total affinity (kcal/mol)
            "minimized_affinity": float(row[0]),
            "rmsd_lb": rmsd_lb,
            "rmsd_ub": rmsd_ub,
            "cnn_score": float("nan"),      # VINA has no CNN scoring
            "cnn_affinity": float("nan"),
        })

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def run_docking(ligand_name: str, config: dict) -> pd.DataFrame:
    """
    Stage 4: Dock one ligand using the configured backend (gnina or vina).

    Reads binding site parameters from data/prepared/binding_site.json.
    Writes per-ligand scores CSV to data/results/{name}_scores.csv.

    Parameters
    ----------
    ligand_name : str
        Must match a row in ligands.csv and have a corresponding .pdbqt file.
    config : dict
        Pipeline config loaded from config.yaml.

    Returns
    -------
    pd.DataFrame
        Scored poses for this ligand.
    """
    site_path = Path(config["paths"]["prepared"]) / "binding_site.json"
    if not site_path.exists():
        raise FileNotFoundError(
            f"Binding site file not found: {site_path}. "
            "Run Stage 3 (define-site) first."
        )
    with open(site_path) as f:
        site = json.load(f)

    logger.info(
        "Binding site: center=(%.3f, %.3f, %.3f), box=(%.1f×%.1f×%.1f) Å, frame=%s",
        site["center_x"], site["center_y"], site["center_z"],
        site["size_x"], site["size_y"], site["size_z"],
        site.get("coordinate_frame", "unknown"),
    )

    backend = config["docking"]["backend"].lower()

    if backend == "gnina":
        df = _run_gnina(ligand_name, site, config)
    elif backend == "gnina_docker":
        df = _run_gnina_docker(ligand_name, site, config)
    elif backend == "vina":
        df = _run_vina(ligand_name, site, config)
    else:
        raise ValueError(
            f"Unknown docking backend: '{backend}'. "
            "Use 'vina', 'gnina' (native, non-macOS), or 'gnina_docker' (Docker, Apple Silicon)."
        )

    # Ensure unified column set
    for col in ["cnn_score", "cnn_affinity", "minimized_affinity", "rmsd_lb", "rmsd_ub"]:
        if col not in df.columns:
            df[col] = float("nan")

    # Save per-ligand scores
    results_dir = Path(config["paths"]["results"])
    scores_path = results_dir / f"{ligand_name}_scores.csv"
    df.to_csv(scores_path, index=False)
    logger.info(
        "Docking complete for %s (%s): %d poses. Scores written to %s",
        ligand_name, backend, len(df), scores_path
    )

    if not df.empty:
        best = df.iloc[0]
        logger.info(
            "  Best pose: score=%.4f, minimized_affinity=%.2f kcal/mol",
            best.get("score", float("nan")),
            best.get("minimized_affinity", float("nan")),
        )

    return df
