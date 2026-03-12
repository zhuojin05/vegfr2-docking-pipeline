"""
Stage 5: Analyse docking results with ProLIF protein-ligand interaction fingerprints.

For each docked ligand:
  1. Converts the top pose PDBQT → PDB (pure Python, no external binary)
  2. Computes ProLIF fingerprint against the cleaned receptor
  3. Identifies key interactions with VEGFR2 ATP-pocket residues:
       Cys919 (hinge H-bond), Glu885, Asp1046 (DFG motif), Phe1047

Plots (data/results/plots/):
  - score_vs_rank.png: score vs pose rank per ligand (both backends)
  - interaction_heatmap.png: ligand × residue interaction heatmap
  - score_comparison.png: best score per ligand, gnina vs vina

Output: data/results/docking_summary.csv
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

# Key VEGFR2 ATP-pocket residues — interactions with these are most
# pharmacologically meaningful for kinase inhibitor validation
KEY_RESIDUES = ["CYS919", "GLU885", "ASP1046", "PHE1047"]


def _pdbqt_to_pdb(pdbqt_path: Path) -> Path | None:
    """
    Convert the first pose in a PDBQT file to PDB format (pure Python).

    PDBQT is a superset of PDB: columns 1-66 are standard PDB fields;
    columns 68+ are AutoDock-specific (partial charge and atom type).
    Truncating those extra columns produces a valid PDB that MDAnalysis
    and ProLIF can read without requiring the obabel external binary.

    Vina output is a multi-model PDBQT (one MODEL/ENDMDL block per pose);
    this function extracts only the first model (best-scoring pose).

    Parameters
    ----------
    pdbqt_path : Path

    Returns
    -------
    Path or None
        PDB path if conversion succeeded, else None.
    """
    pdb_path = pdbqt_path.with_suffix(".pdb")
    try:
        with open(pdbqt_path) as fh_in, open(pdb_path, "w") as fh_out:
            in_first_model = False
            for line in fh_in:
                record = line[:6].strip()
                if record == "MODEL":
                    if in_first_model:
                        break  # second pose — only want first
                    in_first_model = True
                    continue  # skip the MODEL line itself
                elif record == "ENDMDL":
                    break  # end of first pose
                elif record in ("ATOM", "HETATM"):
                    # Truncate PDBQT-specific columns; PDB standard is cols 1–66
                    fh_out.write(line[:66].rstrip() + "\n")
                elif record == "TER":
                    fh_out.write(line[:27].rstrip() + "\n")
        return pdb_path if pdb_path.stat().st_size > 0 else None
    except OSError as exc:
        logger.warning(
            "PDBQT→PDB conversion failed for %s: %s", pdbqt_path.name, exc
        )
        return None


def _compute_prolif(receptor_pdb: Path, pose_pdb: Path, ligand_name: str) -> pd.DataFrame:
    """
    Compute ProLIF protein-ligand interaction fingerprint.

    ProLIF encodes the presence/absence of interaction types per residue as
    a binary fingerprint. Interaction types used here are standard for kinases:
      - HBDonor/HBAcceptor: hydrogen bonds (key for hinge region, e.g. Cys919)
      - Hydrophobic: van der Waals packing in the hydrophobic pocket
      - PiStacking: aromatic stacking (common with Phe1047)
      - PiCation: cation-pi (charged inhibitor moieties)
      - Anionic/Cationic: ionic contacts
      - VdWContact: close contacts not classified above

    Uses RDKit's PDB reader directly rather than MDAnalysis, because MDAnalysis
    lacks explicit bond records for PDB files and falls back to distance-based
    bond guessing, which produces impossible valences for large receptor structures.

    Parameters
    ----------
    receptor_pdb : Path
    pose_pdb : Path
    ligand_name : str

    Returns
    -------
    pd.DataFrame
        Residue-level interaction presence (1/0).
    """
    try:
        from rdkit import Chem
        import prolif
    except ImportError as exc:
        raise ImportError(
            "RDKit and ProLIF required for interaction analysis.\n"
            "Install: pip install prolif  (rdkit is a prolif dependency)"
        ) from exc

    # RDKit's PDB reader uses residue templates for bond assignment — far more
    # reliable than MDAnalysis distance-based guessing for large proteins.
    # sanitize=False defers valence checking; catchErrors in SanitizeMol lets
    # non-standard atoms (e.g. from pdbfixer) pass without crashing.
    rdkit_receptor = Chem.MolFromPDBFile(
        str(receptor_pdb), sanitize=False, removeHs=False
    )
    if rdkit_receptor is None:
        logger.warning("RDKit could not parse receptor PDB: %s", receptor_pdb.name)
        return pd.DataFrame()
    Chem.SanitizeMol(rdkit_receptor, catchErrors=True)
    mol_receptor = prolif.Molecule.from_rdkit(rdkit_receptor)

    rdkit_ligand = Chem.MolFromPDBFile(
        str(pose_pdb), sanitize=False, removeHs=False
    )
    if rdkit_ligand is None:
        logger.warning("RDKit could not parse pose PDB: %s", pose_pdb.name)
        return pd.DataFrame()
    Chem.SanitizeMol(rdkit_ligand, catchErrors=True)
    mol_ligand = prolif.Molecule.from_rdkit(rdkit_ligand)

    fp = prolif.Fingerprint(
        ["HBDonor", "HBAcceptor", "Hydrophobic",
         "PiStacking", "PiCation", "Anionic", "Cationic", "VdWContact"]
    )

    try:
        fp.run_from_iterable([mol_ligand], mol_receptor)
        df = fp.to_dataframe()
        df.index = [ligand_name]
        return df
    except Exception as e:
        logger.warning("ProLIF fingerprint failed for %s: %s", ligand_name, e)
        return pd.DataFrame()


def _plot_score_vs_rank(all_scores: dict[str, pd.DataFrame], plots_dir: Path) -> None:
    """
    Line plot of docking score vs pose rank for each ligand and backend.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    for ligand_name, df in all_scores.items():
        for backend, grp in df.groupby("backend"):
            grp_sorted = grp.sort_values("pose")
            ax.plot(
                grp_sorted["pose"],
                grp_sorted["score"],
                marker="o", label=f"{ligand_name} ({backend})"
            )

    ax.set_xlabel("Pose rank")
    ax.set_ylabel("Score (CNNscore or affinity kcal/mol)")
    ax.set_title("Docking score vs pose rank")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = plots_dir / "score_vs_rank.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    logger.info("Saved %s", out)


def _plot_score_comparison(summary_df: pd.DataFrame, plots_dir: Path) -> None:
    """
    Grouped bar chart: best Gnina vs VINA score per ligand.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    ligands = summary_df["ligand"].tolist()
    x = np.arange(len(ligands))
    width = 0.35

    gnina_scores = summary_df["gnina_best_score"].fillna(0).tolist()
    vina_scores = summary_df["vina_best_score"].fillna(0).tolist()

    ax.bar(x - width / 2, gnina_scores, width, label="Gnina (CNNscore)")
    ax.bar(x + width / 2, vina_scores, width, label="VINA (kcal/mol)")
    ax.set_xticks(x)
    ax.set_xticklabels(ligands, rotation=15, ha="right")
    ax.set_ylabel("Best docking score")
    ax.set_title("Best docking score comparison: Gnina vs VINA")
    ax.legend()
    fig.tight_layout()
    out = plots_dir / "score_comparison.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    logger.info("Saved %s", out)


def _plot_interaction_heatmap(fp_df: pd.DataFrame, plots_dir: Path) -> None:
    """
    Heatmap of ProLIF interaction fingerprints: ligands × residues.
    """
    if fp_df.empty:
        logger.warning("No ProLIF data available — skipping interaction heatmap")
        return

    fig, ax = plt.subplots(figsize=(max(10, fp_df.shape[1] // 2), 4))
    sns.heatmap(
        fp_df.astype(float),
        ax=ax,
        cmap="Blues",
        linewidths=0.5,
        cbar_kws={"label": "Interaction present"},
    )
    ax.set_title("ProLIF interaction fingerprints")
    ax.set_xlabel("Residue:interaction type")
    ax.set_ylabel("Ligand")
    plt.xticks(rotation=45, ha="right", fontsize=6)
    fig.tight_layout()
    out = plots_dir / "interaction_heatmap.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    logger.info("Saved %s", out)


def analyse_results(config: dict) -> pd.DataFrame:
    """
    Stage 5: Compute ProLIF fingerprints and generate summary plots.

    Parameters
    ----------
    config : dict
        Pipeline config loaded from config.yaml.

    Returns
    -------
    pd.DataFrame
        docking_summary.csv contents.
    """
    prepared_dir = Path(config["paths"]["prepared"])
    results_dir = Path(config["paths"]["results"])
    ligands_dir = Path(config["paths"]["ligands"])
    plots_dir = results_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    receptor_pdb = prepared_dir / "4ASD_clean.pdb"
    if not receptor_pdb.exists():
        raise FileNotFoundError(
            f"Cleaned receptor not found: {receptor_pdb}. Run Stage 2 first."
        )

    # Gather all per-ligand score CSVs
    score_files = list(results_dir.glob("*_scores.csv"))
    if not score_files:
        logger.warning("No score CSV files found in %s — run Stage 4 first", results_dir)
        return pd.DataFrame()

    all_scores: dict[str, pd.DataFrame] = {}
    for f in score_files:
        ligand_name = f.stem.replace("_scores", "")
        df = pd.read_csv(f)
        all_scores[ligand_name] = df

    # Collect ProLIF fingerprints
    fp_frames = []
    summary_records = []

    # Load ligand properties for Lipinski flag
    props_path = ligands_dir / "ligand_properties.csv"
    props_df = pd.read_csv(props_path) if props_path.exists() else pd.DataFrame()

    for ligand_name, df in all_scores.items():
        gnina_df = df[df["backend"] == "gnina"]
        vina_df = df[df["backend"] == "vina"]

        gnina_best = gnina_df["score"].max() if not gnina_df.empty else float("nan")
        vina_best = vina_df["score"].min() if not vina_df.empty else float("nan")
        # VINA: more negative = better; Gnina CNNscore: higher = better

        # ProLIF analysis: use top Gnina pose if available, else top VINA pose
        interacting_residues = []
        backend_for_prolif = "gnina" if not gnina_df.empty else "vina"
        pose_suffix = f"{backend_for_prolif}_poses"
        top_pdbqt = results_dir / f"{ligand_name}_{pose_suffix}.pdbqt"

        if top_pdbqt.exists():
            top_pdb = _pdbqt_to_pdb(top_pdbqt)
            if top_pdb:
                fp_df = _compute_prolif(receptor_pdb, top_pdb, ligand_name)
                if not fp_df.empty:
                    fp_frames.append(fp_df)
                    # Extract residues with any interaction
                    # ProLIF columns are MultiIndex: (residue, interaction_type)
                    active_cols = [
                        col for col in fp_df.columns
                        if fp_df[col].iloc[0] == 1 or fp_df[col].iloc[0] is True
                    ]
                    interacting_residues = list(
                        {str(c[0]) if isinstance(c, tuple) else str(c)
                         for c in active_cols}
                    )
                    logger.info(
                        "%s interacting residues: %s",
                        ligand_name, interacting_residues[:10]
                    )

        # Lipinski flag
        lipinski = None
        if not props_df.empty:
            row = props_df[props_df["name"] == ligand_name]
            if not row.empty:
                lipinski = bool(row["Lipinski_pass"].iloc[0])

        summary_records.append({
            "ligand": ligand_name,
            "gnina_best_score": gnina_best,
            "vina_best_score": vina_best,
            "top_interacting_residues": "; ".join(interacting_residues[:10]),
            "lipinski_pass": lipinski,
        })

    summary_df = pd.DataFrame(summary_records)
    summary_path = results_dir / "docking_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    logger.info("Docking summary written to %s", summary_path)

    # Plots
    if all_scores:
        _plot_score_vs_rank(all_scores, plots_dir)
        _plot_score_comparison(summary_df, plots_dir)

    if fp_frames:
        combined_fp = pd.concat(fp_frames)
        _plot_interaction_heatmap(combined_fp, plots_dir)

    return summary_df
