"""Primer splitting: Excel/CSV → V-forward + J-reverse FASTA files."""

import csv
import sys
from pathlib import Path
from typing import Iterator, Tuple


def iter_primer_rows(file_path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (primer_id, sequence) from csv/xlsx first two columns."""
    suffix = file_path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "Reading Excel requires pandas and openpyxl. "
                "Install with: pip install pandas openpyxl"
            ) from exc

        df = pd.read_excel(file_path, engine="openpyxl")
        for _, row in df.iloc[:, :2].iterrows():
            primer_id = str(row.iloc[0]).strip() if len(row) > 0 else ""
            sequence = str(row.iloc[1]).strip() if len(row) > 1 else ""
            yield primer_id, sequence
        return

    if suffix == ".csv":
        with open(file_path, "r", encoding="utf-8-sig", newline="") as f_in:
            reader = csv.reader(f_in)
            for row in reader:
                if not row or len(row) < 2:
                    continue
                yield row[0].strip(), row[1].strip()
        return

    raise ValueError(f"Unsupported file type: {file_path.suffix}")


def split_primers(
    input_path: Path,
    output_dir: Path,
    target: str,
    v_pattern: str = "V",
    j_pattern: str = "J",
) -> Tuple[Path, Path, int, int]:
    """Split primer Excel/CSV into V-forward and J-reverse FASTA files.

    Args:
        input_path: Path to Excel (.xlsx) or CSV file.
        output_dir: Directory for output FASTA files.
        target: Target locus name (e.g., TRB, IGH).
        v_pattern: Substring to identify V-region primers.
        j_pattern: Substring to identify J-region primers.

    Returns:
        Tuple of (v_fasta_path, j_fasta_path, v_count, j_count).
    """
    target_upper = target.upper()
    output_dir.mkdir(parents=True, exist_ok=True)

    v_output = output_dir / f"{target_upper}V_primers.fasta"
    j_output = output_dir / f"{target_upper}J_primers.fasta"

    v_count = 0
    j_count = 0

    with open(v_output, "w", encoding="utf-8") as f_v, \
         open(j_output, "w", encoding="utf-8") as f_j:
        for primer_id, sequence in iter_primer_rows(input_path):
            if not primer_id or not sequence:
                continue
            if primer_id.upper() == "ID" or sequence.upper() == "SEQUENCE":
                continue

            clean_id = primer_id.replace("/", "_").replace(" ", "_")

            if v_pattern.upper() in clean_id.upper():
                f_v.write(f">{clean_id}_fw\n{sequence}\n")
                v_count += 1
            elif j_pattern.upper() in clean_id.upper():
                f_j.write(f">{clean_id}_rev\n{sequence}\n")
                j_count += 1
            else:
                print(f"Warning: cannot determine direction, skipped -> {primer_id}",
                      file=sys.stderr)

    return v_output, j_output, v_count, j_count
