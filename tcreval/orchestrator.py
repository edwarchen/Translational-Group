"""R subprocess orchestration for TCR_Evaluation pipeline."""

import logging
import subprocess
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def _run_rscript(args: list, cwd: Path, timeout: int = 600) -> subprocess.CompletedProcess:
    """Run an R script via Rscript and return the result.

    Args:
        args: Arguments to pass to Rscript (script path + script args).
        cwd: Working directory (project root).
        timeout: Maximum runtime in seconds.

    Returns:
        CompletedProcess with captured stdout/stderr.

    Raises:
        subprocess.TimeoutExpired: If the script takes too long.
        RuntimeError: If the script exits with non-zero code.
    """
    cmd = ["Rscript"] + args
    logger.info("Running: %s", " ".join(str(a) for a in cmd))

    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout,
    )

    if result.stdout:
        for line in result.stdout.strip().split("\n"):
            logger.info("[R] %s", line)

    if result.stderr:
        for line in result.stderr.strip().split("\n"):
            logger.warning("[R] %s", line)

    if result.returncode != 0:
        raise RuntimeError(
            f"R script failed with exit code {result.returncode}:\n{result.stderr}"
        )

    return result


def run_evaluation(
    target: str,
    project_root: Path,
    max_mismatch: int = 3,
    allowed_binding_ratio: float = 1.0,
) -> int:
    """Run the primer coverage evaluation step.

    Args:
        target: Target locus (TRB, IGH, IGK, IGL).
        project_root: Project root directory.
        max_mismatch: Maximum allowed mismatches (default: 3).
        allowed_binding_ratio: Allowed off-target binding ratio (default: 1.0).

    Returns:
        Exit code (0 for success).
    """
    script = project_root / "scripts" / "evaluation.R"
    if not script.exists():
        raise FileNotFoundError(f"Evaluation script not found: {script}")

    logger.info("Running coverage evaluation for target: %s", target)
    _run_rscript([str(script), target], cwd=project_root)
    logger.info("Evaluation complete for %s", target)
    return 0


def run_heatmaps(
    binding_csv: Path,
    ref_fasta: Path,
    output_dir: Path,
    project_root: Path,
    max_plots: Optional[int] = None,
) -> int:
    """Generate per-primer binding heatmaps.

    Args:
        binding_csv: Path to binding sites CSV file.
        ref_fasta: Path to reference FASTA file.
        output_dir: Output directory for PNG files.
        project_root: Project root directory.
        max_plots: Maximum number of plots to generate (default: unlimited).

    Returns:
        Exit code (0 for success).
    """
    script = project_root / "scripts" / "Primer_coverage_heatmap.R"
    if not script.exists():
        raise FileNotFoundError(f"Heatmap script not found: {script}")

    args = [
        str(script),
        f"--binding={binding_csv}",
        f"--ref={ref_fasta}",
        f"--out={output_dir}",
    ]
    if max_plots is not None:
        args.append(f"--max_plots={max_plots}")

    logger.info("Generating primer heatmaps...")
    _run_rscript(args, cwd=project_root)
    logger.info("Heatmaps saved to %s", output_dir)
    return 0


def run_uncovered_plot(
    uncovered_dir: Path,
    project_root: Path,
) -> int:
    """Generate heatmaps for uncovered template sequences.

    Args:
        uncovered_dir: Directory containing uncovered sequence CSV files.
        project_root: Project root directory.

    Returns:
        Exit code (0 for success).
    """
    script = project_root / "scripts" / "Uncovered_plot.R"
    if not script.exists():
        raise FileNotFoundError(f"Uncovered plot script not found: {script}")

    logger.info("Generating uncovered sequence heatmaps...")
    _run_rscript([str(script), f"--dir={uncovered_dir}"], cwd=project_root)
    logger.info("Uncovered plots saved to %s", uncovered_dir)
    return 0


def run_full_pipeline(
    input_file: Path,
    target: str,
    project_root: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    max_mismatch: int = 3,
    allowed_binding_ratio: float = 1.0,
) -> int:
    """Run the complete primer evaluation pipeline end-to-end.

    Args:
        input_file: Path to primer Excel/CSV file.
        target: Target locus (TRB, IGH, IGK, IGL).
        project_root: Project root. Auto-detected if None.
        output_dir: Results directory. Defaults to project_root/results.
        max_mismatch: Max allowed mismatches.
        allowed_binding_ratio: Off-target binding ratio.

    Returns:
        Exit code (0 for success).
    """
    from tcreval.splitter import split_primers
    from tcreval.utils import find_project_root

    if project_root is None:
        project_root = find_project_root()

    if output_dir is None:
        output_dir = project_root / "results"

    target_upper = target.upper()

    # Step 1: Split primers
    logger.info("=== Step 1/4: Splitting primers ===")
    primer_dir = project_root / "data" / "primer_set"
    v_fasta, j_fasta, v_count, j_count = split_primers(
        input_file, primer_dir, target_upper
    )
    logger.info("V primers: %d → %s", v_count, v_fasta)
    logger.info("J primers: %d → %s", j_count, j_fasta)

    # Step 2: Evaluate coverage
    logger.info("=== Step 2/4: Evaluating coverage ===")
    run_evaluation(target_upper, project_root, max_mismatch, allowed_binding_ratio)

    # Step 3: Heatmaps for V and J
    logger.info("=== Step 3/4: Generating primer heatmaps ===")
    cov_dir = output_dir / "coverage_tables"
    heatmap_dir = output_dir / "primer_reference_heatmaps"

    for region in ["V", "J"]:
        binding_csv = cov_dir / f"{target_upper}{region}_binding_sites.csv"
        ref_fasta = project_root / "data" / f"{target_upper}{region}.fasta"
        out_dir = heatmap_dir / f"{target_upper}{region}"

        if binding_csv.exists() and ref_fasta.exists():
            run_heatmaps(binding_csv, ref_fasta, out_dir, project_root)
        else:
            logger.warning("Skipping heatmaps for %s: missing input files", region)

    # Step 4: Uncovered plots
    logger.info("=== Step 4/4: Plotting uncovered sequences ===")
    uncovered_base = output_dir / "uncovered_templates"

    for region in ["V", "J"]:
        unc_dir = uncovered_base / f"{target_upper}{region}_unc"
        if unc_dir.exists() and list(unc_dir.glob("*.csv")):
            run_uncovered_plot(unc_dir, project_root)
        else:
            logger.info("No uncovered sequences for %s%s, skipping plot", target_upper, region)

    logger.info("=== Pipeline complete ===")
    logger.info("Results saved to: %s", output_dir)
    return 0
