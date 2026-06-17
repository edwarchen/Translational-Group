"""CSV data loader and schema validation for TCR report.

Each loader function reads a specific CSV, validates required columns,
and returns a cleaned DataFrame ready for charting.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


# ── Column schemas ──────────────────────────────────────────────────────
# Format: {column_name: (required, default_value)}

SCHEMA_SAMPLE_INFO = {
    "SampleID": (True, None),
    "PatientID": (True, None),
    "Group": (True, None),
    "Type": (False, "Tissue"),
    "SampleType": (False, "Unknown"),
}

SCHEMA_QC_TABLE = {
    "Sample": (True, None),
    "rawYield": (False, 0.0),
    "cleanYield": (False, 0.0),
    "rawReads": (True, None),
    "cleanReads": (True, None),
    "cleanQ30": (True, None),
    "cleanRatio": (True, None),
    "clonalReads": (True, None),
    "clonalRatio": (False, 0.0),
}

SCHEMA_CDR3_LENGTH = {
    "Sample": (True, None),
    "Length": (True, None),
    "OverallFrequency": (True, None),
    "Clone_ID": (True, None),
    "CloneFrequency": (True, None),
}

SCHEMA_V_GENE_FREQ = {
    "Sample": (True, None),
    "V_gene": (True, None),
    "Frequency": (True, None),
}

SCHEMA_J_GENE_FREQ = {
    "Sample": (True, None),
    "J_gene": (True, None),
    "Frequency": (True, None),
}

SCHEMA_VJ_PAIRING = {
    "Sample": (True, None),
    "V_gene": (True, None),
    "J_gene": (True, None),
    "Frequency": (True, None),
}

SCHEMA_CLONE_DIVERSITY = {
    "Sample": (True, None),
    "CloneReads": (True, None),
    "SamplingReads": (False, 1_000_000),
    "CloneCount": (True, None),
    "TopFreq": (True, None),
    "Clonality": (True, None),
    "ShannonIndex": (True, None),
    "SimpsonIndex": (True, None),
    "Evenness": (True, None),
}

SCHEMA_CLONE_COMPOSITION = {
    "Sample": (True, None),
    "FrequencyRange": (True, None),
    "Proportion": (True, None),
}

SCHEMA_FUNCTIONAL_ANNOTATION = {
    "Sample": (True, None),
    "count": (True, None),
    "freq": (True, None),
    "cdr3nt": (False, ""),
    "cdr3aa": (True, None),
    "v": (False, ""),
    "d": (False, ""),
    "j": (False, ""),
    "antigen.species": (True, None),
    "antigen.epitope": (False, ""),
    "antigen.gene": (False, ""),
}

SCHEMA_CANCER_RISK = {
    "sample_id": (True, None),
    "filtered_records": (False, 0),
    "KRAS_type": (False, 0.0),
    "KRAS_score": (False, 0.0),
    "EGFR_type": (False, 0.0),
    "EGFR_score": (False, 0.0),
    "REF_type": (False, 0.0),
    "REF_score": (False, 0.0),
    "LUNG_CANCER_GDNA_type": (False, 0.0),
    "LUNG_CANCER_GDNA_score": (False, 0.0),
    "LUNG_CANCER_TISSUE_type": (False, 0.0),
    "LUNG_CANCER_TISSUE_score": (False, 0.0),
}

SCHEMA_TCR_MOTIF = {
    "Sample": (True, None),
    "Clone_ID": (True, None),
    "CDR3_AA": (True, None),
    "Position": (True, None),
    "AminoAcid": (True, None),
    "Frequency": (True, None),
}

SCHEMA_ANTIGEN_SUMMARY = {
    "Sample": (True, None),
    "antigen": (True, None),
    "freq": (True, None),
    "count": (False, 0),
}

SCHEMA_SIGNIFICANT_EXPANSION = {
    "name": (True, None),
    "sample_Ref": (True, None),
    "sample_Amp": (True, None),
    "fish_pvalue": (True, None),
    "p.adj": (True, None),
    "FC": (True, None),
    "count_Ref": (False, 0),
    "freq_Ref": (True, None),
    "count_Amp": (False, 0),
    "freq_Amp": (True, None),
}


# ── Loader Functions ────────────────────────────────────────────────────

def _validate_and_load(
    filepath: Path,
    schema: dict[str, tuple[bool, any]],
    sample_col: str | None = None,
) -> pd.DataFrame:
    """Load a CSV and validate against schema.

    Parameters
    ----------
    filepath : Path
        Path to the CSV file.
    schema : dict
        Column schema: ``{col: (required, default)}``.
    sample_col : str or None
        If provided, drop rows where this column is NA.

    Returns
    -------
    pd.DataFrame
        Validated DataFrame with required columns.

    Raises
    ------
    FileNotFoundError
        If the CSV file does not exist.
    ValueError
        If a required column is missing.
    """
    if not filepath.exists():
        raise FileNotFoundError(f"Data file not found: {filepath}")

    df = pd.read_csv(filepath)

    # Check required columns
    for col, (required, default) in schema.items():
        if col not in df.columns:
            if required:
                raise ValueError(
                    f"Required column '{col}' missing in {filepath.name}"
                )
            else:
                df[col] = default

    # Keep only schema columns (plus any extras — warn but don't error)
    keep_cols = [c for c in schema if c in df.columns]
    df = df[keep_cols].copy()

    # Drop rows with missing sample identifiers
    if sample_col and sample_col in df.columns:
        df = df.dropna(subset=[sample_col])

    return df


class TCRDataLoader:
    """Loads all TCR data CSVs from a data directory.

    Usage::

        loader = TCRDataLoader("tests/test_data")
        df_sample = loader.load_sample_info()
        df_qc = loader.load_qc_table()
        ...
    """

    def __init__(self, data_dir: str | Path):
        self.data_dir = Path(data_dir)
        if not self.data_dir.is_dir():
            raise NotADirectoryError(f"Data directory not found: {self.data_dir}")

    def _path(self, filename: str) -> Path:
        return self.data_dir / filename

    # ── Individual loaders ──────────────────────────────────────────

    def load_sample_info(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("sample_info.csv"), SCHEMA_SAMPLE_INFO, sample_col="SampleID"
        )

    def load_qc_table(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("qc_table.csv"), SCHEMA_QC_TABLE, sample_col="Sample"
        )

    def load_cdr3_length(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("cdr3_length.csv"), SCHEMA_CDR3_LENGTH, sample_col="Sample"
        )

    def load_v_gene_freq(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("v_gene_freq.csv"), SCHEMA_V_GENE_FREQ, sample_col="Sample"
        )

    def load_j_gene_freq(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("j_gene_freq.csv"), SCHEMA_J_GENE_FREQ, sample_col="Sample"
        )

    def load_vj_pairing(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("vj_pairing.csv"), SCHEMA_VJ_PAIRING, sample_col="Sample"
        )

    def load_clone_diversity(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("clone_diversity.csv"), SCHEMA_CLONE_DIVERSITY, sample_col="Sample"
        )

    def load_clone_composition(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("clone_composition.csv"), SCHEMA_CLONE_COMPOSITION, sample_col="Sample"
        )

    def load_functional_annotation(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("functional_annotation.csv"),
            SCHEMA_FUNCTIONAL_ANNOTATION,
            sample_col="Sample",
        )

    def load_cancer_risk(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("cancer_risk.csv"), SCHEMA_CANCER_RISK, sample_col="sample_id"
        )

    def load_tcr_motif(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("tcr_motif.csv"), SCHEMA_TCR_MOTIF, sample_col="Sample"
        )

    def load_significant_expansion(self) -> pd.DataFrame:
        return _validate_and_load(
            self._path("significant_expansion.csv"),
            SCHEMA_SIGNIFICANT_EXPANSION,
            sample_col="name",
        )

    def load_antigen_summary(self) -> pd.DataFrame:
        path = self._path("antigen_summary.csv")
        if not path.exists():
            return pd.DataFrame(columns=["Sample", "antigen", "freq", "count"])
        return _validate_and_load(path, SCHEMA_ANTIGEN_SUMMARY, sample_col="Sample")

    # ── Convenience ──────────────────────────────────────────────────

    def load_all(self) -> dict[str, pd.DataFrame]:
        """Load all data tables at once.

        Returns a dict mapping table name to DataFrame.
        """
        return {
            "sample_info": self.load_sample_info(),
            "qc_table": self.load_qc_table(),
            "cdr3_length": self.load_cdr3_length(),
            "v_gene_freq": self.load_v_gene_freq(),
            "j_gene_freq": self.load_j_gene_freq(),
            "vj_pairing": self.load_vj_pairing(),
            "clone_diversity": self.load_clone_diversity(),
            "clone_composition": self.load_clone_composition(),
            "functional_annotation": self.load_functional_annotation(),
            "cancer_risk": self.load_cancer_risk(),
            "tcr_motif": self.load_tcr_motif(),
            "significant_expansion": self.load_significant_expansion(),
            "antigen_summary": self.load_antigen_summary(),
        }

    @property
    def sample_ids(self) -> list[str]:
        """Return all sample IDs from the sample_info table."""
        df = self.load_sample_info()
        return df["SampleID"].tolist()
