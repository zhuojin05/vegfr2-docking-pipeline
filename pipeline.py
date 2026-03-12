#!/usr/bin/env python3
"""
VEGFR2 Docking Pipeline — CLI orchestrator.

Usage examples:
  python pipeline.py prepare-ligands --input data/ligands/ligands.csv
  python pipeline.py prepare-receptor
  python pipeline.py define-site
  python pipeline.py dock --ligand sorafenib
  python pipeline.py dock --all
  python pipeline.py analyse
  python pipeline.py run-all --input data/ligands/ligands.csv
"""

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

import yaml


def _setup_logging(log_level: str, log_dir: Path) -> None:
    """Configure logging to console and timestamped file."""
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"pipeline_{timestamp}.log"

    level = getattr(logging, log_level.upper(), logging.INFO)
    fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"

    logging.basicConfig(
        level=level,
        format=fmt,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file),
        ],
    )
    logging.getLogger(__name__).info("Log file: %s", log_file)


def _load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Sub-command handlers
# ---------------------------------------------------------------------------

def cmd_prepare_ligands(args, config: dict) -> None:
    from pipeline.prepare_ligands import prepare_ligands

    csv_path = Path(args.input)
    if not csv_path.exists():
        logging.error("Ligand CSV not found: %s", csv_path)
        sys.exit(1)

    df = prepare_ligands(csv_path, config)
    if df.empty:
        logging.error("No ligands were successfully processed")
        sys.exit(1)
    print(df.to_string(index=False))


def cmd_prepare_receptor(args, config: dict) -> None:
    from pipeline.prepare_receptor import prepare_receptor

    clean_pdb = prepare_receptor(config)
    logging.info("Receptor preparation complete: %s", clean_pdb)


def cmd_define_site(args, config: dict) -> None:
    from pipeline.define_binding_site import define_binding_site

    site = define_binding_site(config)
    print("\nBinding site parameters:")
    for k, v in site.items():
        print(f"  {k}: {v}")


def cmd_dock(args, config: dict) -> None:
    from pipeline.run_docking import run_docking
    import pandas as pd

    if args.all:
        # Load all ligand names from the ligands CSV
        ligands_csv = Path(config["paths"]["ligands"]) / "ligands.csv"
        if not ligands_csv.exists():
            # Try positional --input if provided
            logging.error("ligands.csv not found at %s", ligands_csv)
            sys.exit(1)
        import pandas as _pd
        names = _pd.read_csv(ligands_csv)["name"].tolist()
    elif args.ligand:
        names = [args.ligand]
    else:
        logging.error("Specify --ligand <name> or --all")
        sys.exit(1)

    all_dfs = []
    for name in names:
        try:
            df = run_docking(name, config)
            all_dfs.append(df)
        except Exception as e:
            logging.warning("Docking failed for %s: %s", name, e)

    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        print(combined.to_string(index=False))


def cmd_analyse(args, config: dict) -> None:
    from pipeline.analyse_results import analyse_results

    summary = analyse_results(config)
    if not summary.empty:
        print(summary.to_string(index=False))


def cmd_run_all(args, config: dict) -> None:
    """Run all pipeline stages sequentially."""
    from pipeline.prepare_ligands import prepare_ligands
    from pipeline.prepare_receptor import prepare_receptor
    from pipeline.define_binding_site import define_binding_site
    from pipeline.run_docking import run_docking
    from pipeline.analyse_results import analyse_results
    import pandas as pd

    log = logging.getLogger("run-all")

    # Stage 1
    log.info("=== Stage 1: Prepare ligands ===")
    csv_path = Path(args.input)
    props_df = prepare_ligands(csv_path, config)
    ligand_names = props_df["name"].tolist() if not props_df.empty else []
    log.info("Processed %d ligands", len(ligand_names))

    # Stage 2
    log.info("=== Stage 2: Prepare receptor ===")
    prepare_receptor(config)

    # Stage 3
    log.info("=== Stage 3: Define binding site ===")
    define_binding_site(config)

    # Stage 4
    log.info("=== Stage 4: Run docking ===")
    for name in ligand_names:
        try:
            run_docking(name, config)
        except Exception as e:
            log.warning("Docking failed for %s — continuing: %s", name, e)

    # Stage 5
    log.info("=== Stage 5: Analyse results ===")
    summary = analyse_results(config)
    if not summary.empty:
        print(summary.to_string(index=False))

    log.info("Pipeline complete.")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="VEGFR2 molecular docking pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config", default="config.yaml",
        help="Path to config.yaml (default: config.yaml)"
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # prepare-ligands
    p_lig = subparsers.add_parser("prepare-ligands", help="Stage 1: SMILES → pdbqt")
    p_lig.add_argument(
        "--input", default="data/ligands/ligands.csv",
        help="Path to ligands CSV (default: data/ligands/ligands.csv)"
    )

    # prepare-receptor
    subparsers.add_parser("prepare-receptor", help="Stage 2: PDB → pdbqt")

    # define-site
    subparsers.add_parser("define-site", help="Stage 3: Extract binding site from 3WZE")

    # dock
    p_dock = subparsers.add_parser("dock", help="Stage 4: Run docking")
    p_dock.add_argument("--ligand", help="Name of one ligand to dock")
    p_dock.add_argument("--all", action="store_true", help="Dock all ligands in ligands.csv")

    # analyse
    subparsers.add_parser("analyse", help="Stage 5: ProLIF analysis + plots")

    # run-all
    p_all = subparsers.add_parser("run-all", help="Run all stages sequentially")
    p_all.add_argument(
        "--input", default="data/ligands/ligands.csv",
        help="Path to ligands CSV (default: data/ligands/ligands.csv)"
    )

    args = parser.parse_args()
    config = _load_config(args.config)
    _setup_logging(args.log_level, Path(config["paths"]["logs"]))

    dispatch = {
        "prepare-ligands": cmd_prepare_ligands,
        "prepare-receptor": cmd_prepare_receptor,
        "define-site": cmd_define_site,
        "dock": cmd_dock,
        "analyse": cmd_analyse,
        "run-all": cmd_run_all,
    }

    dispatch[args.command](args, config)


if __name__ == "__main__":
    main()
