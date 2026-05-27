"""Path resolution and utilities for TCR_Evaluation."""

import logging
import sys
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    """Configure logging based on verbosity flags."""
    if quiet:
        level = logging.WARNING
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


def resolve_path(path: str, root: Optional[Path] = None) -> Path:
    """Resolve a path relative to project root or CWD.

    Args:
        path: File path (absolute or relative).
        root: Project root directory. If None, uses CWD.

    Returns:
        Resolved absolute Path.
    """
    p = Path(path)
    if p.is_absolute():
        return p

    if root is None:
        root = Path.cwd()

    return (root / p).resolve()


def check_r_available() -> bool:
    """Verify Rscript is on PATH and openPrimeR is installed."""
    import subprocess

    try:
        result = subprocess.run(
            ["Rscript", "-e", "library(openPrimeR)"],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            logger.error(
                "openPrimeR is not installed in R.\n"
                "Install it with:\n"
                '  R -e \'if (!require("BiocManager")) install.packages("BiocManager"); '
                'BiocManager::install(c("openPrimeR", "Biostrings"))\''
            )
            return False
        return True
    except FileNotFoundError:
        logger.error("Rscript not found. Please install R 4.0+ first.")
        return False


def find_project_root() -> Path:
    """Locate the project root by searching for R/core.R or config/."""
    current = Path.cwd().resolve()
    for parent in [current] + list(current.parents):
        if (parent / "R" / "core.R").exists():
            return parent
        if (parent / "config" / "defaults.yaml").exists():
            return parent
    return current
