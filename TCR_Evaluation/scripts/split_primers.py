#!/usr/bin/env python3
"""Split a primer Excel/CSV file into V-forward and J-reverse FASTA files."""

import argparse
import csv
import sys
from pathlib import Path


def iter_primer_rows(file_path: Path):
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


def main():
    parser = argparse.ArgumentParser(
        description="Split primer Excel/CSV into V-forward and J-reverse FASTA files."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to input Excel (.xlsx) or CSV file with columns: Primer ID, Sequence"
    )
    parser.add_argument(
        "-t", "--target", required=True,
        help="Target locus name (e.g., TRB, IGH, IGK, IGL). Used to name output files."
    )
    parser.add_argument(
        "-o", "--output-dir", default="data/primer_set",
        help="Output directory for FASTA files (default: data/primer_set)"
    )
    parser.add_argument(
        "--v-pattern", default="V",
        help="Substring to identify V-region primers in ID (default: V)"
    )
    parser.add_argument(
        "--j-pattern", default="J",
        help="Substring to identify J-region primers in ID (default: J)"
    )
    args = parser.parse_args()

    src = Path(args.input)
    if not src.exists():
        print(f"Error: input file not found: '{args.input}'", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    target_upper = args.target.upper()
    v_output = out_dir / f"{target_upper}V_primers.fasta"
    j_output = out_dir / f"{target_upper}J_primers.fasta"

    v_count = 0
    j_count = 0

    with open(v_output, "w", encoding="utf-8") as f_v, \
         open(j_output, "w", encoding="utf-8") as f_j:
        for primer_id, sequence in iter_primer_rows(src):
            if not primer_id or not sequence:
                continue
            if primer_id.upper() == "ID" or sequence.upper() == "SEQUENCE":
                continue

            clean_id = primer_id.replace("/", "_").replace(" ", "_")

            if args.v_pattern.upper() in clean_id.upper():
                f_v.write(f">{clean_id}_fw\n{sequence}\n")
                v_count += 1
            elif args.j_pattern.upper() in clean_id.upper():
                f_j.write(f">{clean_id}_rev\n{sequence}\n")
                j_count += 1
            else:
                print(f"Warning: cannot determine direction for primer, skipped -> {primer_id}")

    print(f"Done. Extracted {v_count} V-forward primers -> {v_output}")
    print(f"       Extracted {j_count} J-reverse primers -> {j_output}")

    if v_count == 0 and j_count == 0:
        print("Warning: No primers extracted. Check your input file and V/J patterns.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
