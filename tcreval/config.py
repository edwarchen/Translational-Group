"""YAML configuration loader for TCR_Evaluation."""

from pathlib import Path
from typing import Any, Dict, Optional

import yaml


def find_project_root() -> Path:
    """Locate the project root by searching for config/defaults.yaml."""
    current = Path.cwd().resolve()
    for parent in [current] + list(current.parents):
        if (parent / "config" / "defaults.yaml").exists():
            return parent
        if (parent / "R" / "core.R").exists():
            return parent
    return current


def load_config(config_path: Optional[Path] = None) -> Dict[str, Any]:
    """Load configuration from YAML file.

    Args:
        config_path: Path to YAML config. If None, uses config/defaults.yaml
                     relative to project root.

    Returns:
        Dict with evaluation settings and target configurations.
    """
    if config_path is None:
        config_path = find_project_root() / "config" / "defaults.yaml"

    if not config_path.exists():
        raise FileNotFoundError(
            f"Config file not found: {config_path}\n"
            "Run from the project root or provide --config."
        )

    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    return config


def get_target_config(target: str, config: Optional[Dict] = None) -> Dict[str, str]:
    """Get paths and settings for a specific target locus.

    Args:
        target: Locus name (TRB, IGH, IGK, IGL).
        config: Full configuration dict. Loaded from defaults.yaml if None.

    Returns:
        Dict with v_ref, j_ref, v_primers, j_primers paths.
    """
    if config is None:
        config = load_config()

    target = target.upper()
    targets = config.get("targets", {})

    if target not in targets:
        available = ", ".join(targets.keys())
        raise ValueError(
            f"Unknown target: {target}. Available targets: {available}"
        )

    return targets[target]


def get_eval_settings(config: Optional[Dict] = None) -> Dict[str, Any]:
    """Get evaluation parameters from config.

    Returns:
        Dict with max_mismatch, allowed_binding_ratio, trim settings.
    """
    if config is None:
        config = load_config()

    return config.get("evaluation", {})
