"""Streamlit web UI for TCR/BCR primer evaluation pipeline.

Usage:
    streamlit run app.py
"""

import io
import shutil
import zipfile
from pathlib import Path

import streamlit as st

# Page config must be the first Streamlit command
st.set_page_config(
    page_title="TCR/BCR Primer Evaluation",
    page_icon="🧬",
    layout="wide",
)


def check_environment() -> dict:
    """Verify R and Python dependencies are available."""
    import subprocess

    status = {"r_ok": False, "r_message": ""}

    try:
        result = subprocess.run(
            ["Rscript", "-e", "library(openPrimeR); library(Biostrings); library(ggplot2); library(tidyr)"],
            capture_output=True, text=True, timeout=30,
        )
        if result.returncode == 0:
            status["r_ok"] = True
            status["r_message"] = "R + openPrimeR ready"
        else:
            status["r_message"] = "openPrimeR not found. Run: Rscript R/install_deps.R"
    except FileNotFoundError:
        status["r_message"] = "Rscript not found. Please install R 4.0+."

    return status


def find_project_root() -> Path:
    """Locate the project root (where R/core.R lives)."""
    return Path(__file__).resolve().parent


def run_pipeline_from_table(
    input_path: Path,
    target: str,
    project_root: Path,
    max_mismatch: int = 3,
    binding_ratio: float = 1.0,
) -> Path:
    """Run full pipeline from an Excel/CSV primer table (split + eval + heatmaps + uncovered)."""
    from tcreval.orchestrator import run_full_pipeline

    results_dir = project_root / "results"
    run_full_pipeline(
        input_file=input_path,
        target=target,
        project_root=project_root,
        output_dir=results_dir,
        max_mismatch=max_mismatch,
        allowed_binding_ratio=binding_ratio,
    )
    return results_dir


def run_evaluation_only(
    target: str,
    project_root: Path,
    max_mismatch: int = 3,
    binding_ratio: float = 1.0,
) -> Path:
    """Run only the evaluation + heatmap + uncovered steps (primers already in FASTA format)."""
    from tcreval.orchestrator import run_evaluation, run_heatmaps, run_uncovered_plot
    import logging

    logging.basicConfig(level=logging.INFO)

    results_dir = project_root / "results"

    # Step 1: Evaluate
    run_evaluation(target, project_root, max_mismatch, binding_ratio)

    # Step 2: Heatmaps
    cov_dir = results_dir / "coverage_tables"
    heatmap_dir = results_dir / "primer_reference_heatmaps"
    for region in ["V", "J"]:
        binding_csv = cov_dir / f"{target}{region}_binding_sites.csv"
        ref_fasta = project_root / "data" / f"{target}{region}.fasta"
        out_dir = heatmap_dir / f"{target}{region}"
        if binding_csv.exists() and ref_fasta.exists():
            run_heatmaps(binding_csv, ref_fasta, out_dir, project_root)

    # Step 3: Uncovered plots
    uncovered_base = results_dir / "uncovered_templates"
    for region in ["V", "J"]:
        unc_dir = uncovered_base / f"{target}{region}_unc"
        if unc_dir.exists() and list(unc_dir.glob("*.csv")):
            run_uncovered_plot(unc_dir, project_root)

    return results_dir


def create_results_zip(results_dir: Path) -> bytes:
    """Create a zip archive of all results."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(results_dir.rglob("*")):
            if f.is_file():
                zf.write(f, f.relative_to(results_dir))
    buf.seek(0)
    return buf.getvalue()


def display_results(results_dir: Path, target: str):
    """Render the four result tabs."""
    st.markdown("---")
    st.header("📊 Results")

    tab_cov, tab_heatmaps, tab_uncovered, tab_download = st.tabs([
        "Coverage", "Primer Heatmaps", "Uncovered Templates", "Download"
    ])

    cov_dir = results_dir / "coverage_tables"
    plot_dir = results_dir / "coverage_plots"
    heatmap_dir = results_dir / "primer_reference_heatmaps"
    uncovered_dir = results_dir / "uncovered_templates"

    with tab_cov:
        col1, col2 = st.columns(2)
        for region in ["V", "J"]:
            region_key = f"{target}{region}"
            cov_plot = plot_dir / f"{region_key}_cov.png"
            cov_csv = cov_dir / f"{region_key}_cov.csv"

            with col1 if region == "V" else col2:
                st.subheader(f"{region}-region Coverage")
                if cov_csv.exists():
                    import pandas as pd
                    df = pd.read_csv(cov_csv)
                    st.dataframe(df, use_container_width=True, hide_index=True)
                if cov_plot.exists():
                    st.image(str(cov_plot), use_container_width=True)

    with tab_heatmaps:
        st.subheader("Per-Primer Binding Alignments")
        for region_tag in [f"{target}V", f"{target}J"]:
            hm_dir = heatmap_dir / region_tag
            if hm_dir.exists():
                pngs = sorted(hm_dir.glob("*.png"))
                if pngs:
                    st.caption(f"{region_tag}: {len(pngs)} primers")
                    cols = st.columns(3)
                    for i, png in enumerate(pngs):
                        with cols[i % 3]:
                            st.image(str(png), caption=png.stem, use_container_width=True)

    with tab_uncovered:
        st.subheader("Uncovered Templates by Gene Family")
        for region_tag in [f"{target}V_unc", f"{target}J_unc"]:
            unc_dir = uncovered_dir / region_tag
            if unc_dir.exists():
                pngs = sorted(unc_dir.glob("*_heatmap.png"))
                if pngs:
                    st.caption(f"{region_tag}: {len(pngs)} families with gaps")
                    cols = st.columns(2)
                    for i, png in enumerate(pngs):
                        with cols[i % 2]:
                            st.image(str(png), caption=png.stem, use_container_width=True)
                else:
                    region_label = region_tag.replace("_unc", "")
                    st.info(f"✅ {region_label}: all families fully covered")

    with tab_download:
        st.subheader("Download Results Package")
        zip_data = create_results_zip(results_dir)
        st.download_button(
            label=f"📥 Download {target}_results.zip",
            data=zip_data,
            file_name=f"{target}_evaluation_results.zip",
            mime="application/zip",
            use_container_width=True,
        )
        st.caption(f"Results directory: `{results_dir}`")


# ═══ Sidebar ═══

st.sidebar.title("🧬 TCR/BCR Primer Evaluation")
st.sidebar.markdown("---")

env = check_environment()
if env["r_ok"]:
    st.sidebar.success("✅ Environment ready")
else:
    st.sidebar.error(f"⚠️ {env['r_message']}")

st.sidebar.markdown("---")

# Input mode
input_mode = st.sidebar.radio(
    "Input Mode",
    options=["Primer Table (.xlsx/.csv)", "Primer FASTA (.fasta/.fa)"],
    help=(
        "**Primer Table**: Excel or CSV with two columns (ID, Sequence). "
        "The tool will automatically split into V-forward and J-reverse FASTA files.\n\n"
        "**Primer FASTA**: Pre-split FASTA files. IDs must end with `_fw` (V-forward) "
        "or `_rev` (J-reverse). The split step is skipped."
    ),
)

# Target selection
target = st.sidebar.selectbox(
    "Target Locus",
    options=["TRB", "IGH", "IGK", "IGL"],
    help="Immune receptor locus to evaluate",
)

# Advanced settings
with st.sidebar.expander("Advanced Settings"):
    max_mismatch = st.slider(
        "Max Mismatches", min_value=0, max_value=5, value=3,
        help="Maximum allowed mismatches for primer binding",
    )
    binding_ratio = st.number_input(
        "Off-target Binding Ratio", min_value=0.0, max_value=2.0, value=1.0, step=0.1,
        help="Allowed ratio of off-target to on-target binding",
    )

st.sidebar.markdown("---")

project_root = find_project_root()

# File upload — depends on mode
if input_mode == "Primer Table (.xlsx/.csv)":
    uploaded_file = st.sidebar.file_uploader(
        "Upload Primer Table",
        type=["xlsx", "xls", "csv"],
        help="Columns: Primer ID, Sequence. IDs containing 'V' → forward, 'J' → reverse.",
    )
else:
    st.sidebar.caption("Upload one combined FASTA, or two separate files (V + J).")
    v_fasta = st.sidebar.file_uploader(
        "V-forward Primers (.fasta)", type=["fasta", "fa", "fna"],
        help="FASTA with IDs ending in _fw (e.g., >TRBV1_fw)",
    )
    j_fasta = st.sidebar.file_uploader(
        "J-reverse Primers (.fasta)", type=["fasta", "fa", "fna"],
        help="FASTA with IDs ending in _rev (e.g., >TRBJ1_rev)",
    )
    combined_fasta = st.sidebar.file_uploader(
        "Or: Combined FASTA (V+J)", type=["fasta", "fa", "fna"],
        help="Single FASTA containing both _fw and _rev primers",
    )

# ═══ Main area ═══

st.title("TCR/BCR Primer Coverage Evaluation")

# --- Primer Table mode ---
if input_mode == "Primer Table (.xlsx/.csv)":
    if not uploaded_file:
        st.info("👈 Upload a primer table in the sidebar to get started.")
        st.markdown("""
        ### What this tool does

        1. **Parses** your primer list into V-forward and J-reverse FASTA files
        2. **Evaluates** primer coverage against IMGT reference sequences using in silico PCR
        3. **Visualizes** binding sites, coverage statistics, and uncovered templates

        ### Input Format
        **Primer Table** (Excel/CSV) — two columns:
        - **Primer ID**: e.g., `TRBV1`, `IGHJ4` (must contain "V" or "J")
        - **Sequence**: e.g., `GCTACTTCGGAGCCTCGG`

        **Primer FASTA** — IDs must end with `_fw` or `_rev`:
        ```
        >TRBV1_fw
        GCTACTTCGGAGCCTCGG
        >TRBJ1_rev
        CTTACCTGAGGAGACGGTGACC
        ```
        """)
        st.stop()

    saved_path = project_root / "data" / "primer_set" / uploaded_file.name
    saved_path.parent.mkdir(parents=True, exist_ok=True)
    with open(saved_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    st.success(f"📄 Uploaded: **{uploaded_file.name}** ({uploaded_file.size:,} bytes)")

    if st.button("🚀 Run Evaluation", type="primary", use_container_width=True):
        if not env["r_ok"]:
            st.error(f"Cannot run: {env['r_message']}")
            st.stop()

        try:
            with st.status("Running pipeline...", expanded=True) as status:
                st.write("**Step 1/4:** Splitting primers into V/J FASTA files...")
                st.write("**Step 2/4:** Evaluating coverage with openPrimeR...")
                st.write("**Step 3/4:** Generating primer binding heatmaps...")
                st.write("**Step 4/4:** Plotting uncovered templates...")

                results_dir = run_pipeline_from_table(
                    input_path=saved_path,
                    target=target,
                    project_root=project_root,
                    max_mismatch=max_mismatch,
                    binding_ratio=binding_ratio,
                )
                status.update(label="Pipeline complete!", state="complete")

            st.success(f"✅ Evaluation complete for **{target}**")
            display_results(results_dir, target)

        except Exception as e:
            st.error(f"Pipeline failed: {e}")

# --- FASTA mode ---
else:
    # Determine which FASTA files to use
    fasta_inputs = []
    if combined_fasta:
        fasta_inputs.append(("combined", combined_fasta))
    if v_fasta:
        fasta_inputs.append(("V", v_fasta))
    if j_fasta:
        fasta_inputs.append(("J", j_fasta))

    if not fasta_inputs:
        st.info("👈 Upload primer FASTA files in the sidebar to get started.")
        st.markdown("""
        ### FASTA Mode

        Upload primers already in FASTA format. Two options:

        1. **Combined FASTA** — one file with both V and J primers (IDs with `_fw` / `_rev` suffixes)
        2. **Separate V/J FASTA** — two files, one for V-forward, one for J-reverse

        The tool will place the FASTA files in the correct location and run evaluation directly,
        skipping the Excel→FASTA split step.
        """)
        st.stop()

    primer_dir = project_root / "data" / "primer_set"
    primer_dir.mkdir(parents=True, exist_ok=True)
    saved_files = []

    for label, f in fasta_inputs:
        if label == "combined":
            # Parse combined FASTA: extract _fw → V file, _rev → J file
            content = f.getvalue().decode("utf-8")
            v_lines = []
            j_lines = []
            current_id = None
            current_lines = []
            target_list = None

            for line in content.split("\n"):
                stripped = line.strip()
                if stripped.startswith(">"):
                    # flush previous entry
                    if current_id and target_list is not None:
                        target_list.append(f">{current_id}\n" + "".join(current_lines))
                    current_id = stripped[1:]
                    current_lines = []
                    if current_id.endswith("_fw"):
                        target_list = v_lines
                    elif current_id.endswith("_rev"):
                        target_list = j_lines
                    else:
                        target_list = None
                elif target_list is not None and stripped:
                    current_lines.append(stripped + "\n")

            # flush last entry
            if current_id and target_list is not None:
                target_list.append(f">{current_id}\n" + "".join(current_lines))

            v_path = primer_dir / f"{target}V_primers.fasta"
            j_path = primer_dir / f"{target}J_primers.fasta"
            with open(v_path, "w") as fv:
                fv.writelines(v_lines)
            with open(j_path, "w") as fj:
                fj.writelines(j_lines)
            saved_files.extend([v_path, j_path])
            st.success(f"📄 Parsed combined FASTA: {len(v_lines)} V + {len(j_lines)} J primers")

        elif label == "V":
            dest = primer_dir / f"{target}V_primers.fasta"
            with open(dest, "wb") as fdest:
                fdest.write(f.getvalue())
            saved_files.append(dest)
            st.success(f"📄 V-forward primers saved ({f.size:,} bytes)")

        elif label == "J":
            dest = primer_dir / f"{target}J_primers.fasta"
            with open(dest, "wb") as fdest:
                fdest.write(f.getvalue())
            saved_files.append(dest)
            st.success(f"📄 J-reverse primers saved ({f.size:,} bytes)")

    if st.button("🚀 Run Evaluation", type="primary", use_container_width=True):
        if not env["r_ok"]:
            st.error(f"Cannot run: {env['r_message']}")
            st.stop()

        try:
            with st.status("Running pipeline...", expanded=True) as status:
                st.write("**Step 1/3:** Evaluating coverage with openPrimeR...")
                st.write("**Step 2/3:** Generating primer binding heatmaps...")
                st.write("**Step 3/3:** Plotting uncovered templates...")

                results_dir = run_evaluation_only(
                    target=target,
                    project_root=project_root,
                    max_mismatch=max_mismatch,
                    binding_ratio=binding_ratio,
                )
                status.update(label="Pipeline complete!", state="complete")

            st.success(f"✅ Evaluation complete for **{target}**")
            display_results(results_dir, target)

        except Exception as e:
            st.error(f"Pipeline failed: {e}")
