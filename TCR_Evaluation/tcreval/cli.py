"""Command-line interface for TCR/BCR primer evaluation.

Usage:
    tcreval run   -i primers.xlsx -t TRB       # Full pipeline
    tcreval split -i primers.xlsx -t TRB       # Split primers only
    tcreval eval  -t IGH                       # Evaluate coverage only
    tcreval heatmap -b bindings.csv -r ref.fasta -o out/
    tcreval uncovered -d results/uncovered_templates/IGHV_unc
"""

import sys
from pathlib import Path

import click

from tcreval import __version__


@click.group()
@click.version_option(version=__version__)
def main():
    """TCR_Evaluation — TCR/BCR multiplex PCR primer coverage evaluation toolkit."""
    pass


@main.command()
@click.option("-i", "--input", required=True, type=click.Path(exists=True),
              help="Input Excel (.xlsx) or CSV file with primer ID and sequence columns.")
@click.option("-t", "--target", required=True,
              help="Target locus: TRB, IGH, IGK, or IGL.")
@click.option("-o", "--output-dir", default="data/primer_set",
              help="Output directory for FASTA files (default: data/primer_set).")
@click.option("--v-pattern", default="V",
              help="Substring to identify V-region primers (default: V).")
@click.option("--j-pattern", default="J",
              help="Substring to identify J-region primers (default: J).")
def split(input, target, output_dir, v_pattern, j_pattern):
    """Split primer Excel/CSV into V-forward and J-reverse FASTA files."""
    from tcreval.splitter import split_primers

    v_out, j_out, v_count, j_count = split_primers(
        Path(input), Path(output_dir), target, v_pattern, j_pattern
    )
    click.echo(f"V primers: {v_count:>3d} → {v_out}")
    click.echo(f"J primers: {j_count:>3d} → {j_out}")


@main.command()
@click.option("-t", "--target", required=True,
              help="Target locus: TRB, IGH, IGK, or IGL.")
@click.option("--max-mismatch", default=3, type=int,
              help="Maximum allowed mismatches (default: 3).")
@click.option("--binding-ratio", default=1.0, type=float,
              help="Allowed off-target binding ratio (default: 1.0).")
def eval(target, max_mismatch, binding_ratio):
    """Run primer coverage evaluation for a target locus."""
    from tcreval.orchestrator import run_evaluation
    from tcreval.utils import find_project_root

    try:
        run_evaluation(target, find_project_root(), max_mismatch, binding_ratio)
    except RuntimeError as e:
        click.echo(str(e), err=True)
        sys.exit(1)


@main.command()
@click.option("-b", "--binding", required=True, type=click.Path(exists=True),
              help="Path to binding sites CSV file.")
@click.option("-r", "--ref", required=True, type=click.Path(exists=True),
              help="Path to reference FASTA file.")
@click.option("-o", "--out", required=True, type=click.Path(),
              help="Output directory for heatmap PNG files.")
@click.option("--max-plots", default=None, type=int,
              help="Maximum number of plots to generate.")
def heatmap(binding, ref, out, max_plots):
    """Generate per-primer nucleotide alignment heatmaps."""
    from tcreval.orchestrator import run_heatmaps
    from tcreval.utils import find_project_root

    try:
        run_heatmaps(Path(binding), Path(ref), Path(out),
                     find_project_root(), max_plots)
    except RuntimeError as e:
        click.echo(str(e), err=True)
        sys.exit(1)


@main.command()
@click.option("-d", "--dir", required=True, type=click.Path(exists=True),
              help="Directory containing uncovered template CSV files.")
def uncovered(dir):
    """Generate heatmaps for uncovered template sequences."""
    from tcreval.orchestrator import run_uncovered_plot
    from tcreval.utils import find_project_root

    try:
        run_uncovered_plot(Path(dir), find_project_root())
    except RuntimeError as e:
        click.echo(str(e), err=True)
        sys.exit(1)


@main.command()
@click.option("-i", "--input", required=True, type=click.Path(exists=True),
              help="Input Excel (.xlsx) or CSV file with primer ID and sequence columns.")
@click.option("-t", "--target", required=True,
              help="Target locus: TRB, IGH, IGK, or IGL.")
@click.option("-o", "--output-dir", default=None, type=click.Path(),
              help="Results directory (default: results/ in project root).")
@click.option("--max-mismatch", default=3, type=int,
              help="Maximum allowed mismatches (default: 3).")
@click.option("--binding-ratio", default=1.0, type=float,
              help="Allowed off-target binding ratio (default: 1.0).")
@click.option("-v", "--verbose", is_flag=True, help="Verbose logging.")
@click.option("-q", "--quiet", is_flag=True, help="Suppress non-error output.")
def run(input, target, output_dir, max_mismatch, binding_ratio, verbose, quiet):
    """Run the complete pipeline: split → evaluate → heatmaps → uncovered plots."""
    from tcreval.orchestrator import run_full_pipeline
    from tcreval.utils import setup_logging

    setup_logging(verbose=verbose, quiet=quiet)

    try:
        exit_code = run_full_pipeline(
            input_file=Path(input),
            target=target,
            output_dir=Path(output_dir) if output_dir else None,
            max_mismatch=max_mismatch,
            allowed_binding_ratio=binding_ratio,
        )
        sys.exit(exit_code)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
