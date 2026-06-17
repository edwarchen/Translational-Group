#!/usr/bin/env python3
"""Generate TCR report CSV tables from the cloud-platform TSV results file.

Combines the aggregate TSV (QC, diversity, cancer risk, antigen annotation,
gene frequencies, AA composition) with per-clone detail from a clonotype
TXT file (CDR3 sequences, V-J pairing, motif data).

Usage:
    python convert_from_tsv.py Results_vdj_tcrmatch.tsv sample.clonotypes.TRB.txt -o output/
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")
RANGE_ORDER = [
    "Rare (≤0.01%)", "Small (0.01%-0.1%)", "Medium (0.1%-1%)",
    "Large (1%-10%)", "Hyperexpanded (>10%)",
]


def clean_gene(raw: str) -> str | None:
    if pd.isna(raw) or not str(raw).strip():
        return None
    first = str(raw).split(",")[0].strip()
    m = re.match(r"(TRB[VDJ]\d+-\d+)", first)
    return m.group(1) if m else None


def compute_diversity(freqs: np.ndarray) -> dict:
    n = len(freqs)
    if n == 0:
        return {"ShannonIndex": 0, "SimpsonIndex": 0, "Evenness": 0, "Clonality": 0}
    shannon = -np.sum(freqs * np.log(freqs + 1e-300))
    simpson = np.sum(freqs ** 2)
    max_shannon = np.log(n) if n > 1 else 1.0
    evenness = shannon / max_shannon if max_shannon > 0 else 0
    return {
        "ShannonIndex": round(float(shannon), 6),
        "SimpsonIndex": round(float(simpson), 6),
        "Evenness": round(float(evenness), 6),
        "Clonality": round(float(1 - evenness), 6),
    }


def classify_freq_range(freq: float) -> str:
    if freq > 0.10:      return "Hyperexpanded (>10%)"
    elif freq > 0.01:    return "Large (1%-10%)"
    elif freq > 0.001:   return "Medium (0.1%-1%)"
    elif freq > 0.0001:  return "Small (0.01%-0.1%)"
    else:                return "Rare (≤0.01%)"


def convert(tsv_path: str, clonotype_path: str, output_dir: str):
    tsv_path = Path(tsv_path)
    clonotype_path = Path(clonotype_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load TSV aggregate data ──────────────────────────────────────────
    print(f"Loading {tsv_path} ...")
    tsv = pd.read_csv(tsv_path, sep="\t")
    print(f"  {len(tsv)} samples found")

    # Load clonotype detail
    print(f"Loading {clonotype_path} ...")
    cl = pd.read_csv(clonotype_path, sep="\t")
    # Filter non-productive
    n_before = len(cl)
    cl = cl[~cl["aaFeatureSequences(CDR3)"].str.contains(r"\*", na=False)]
    cl = cl[~cl["aaFeatureSequences(CDR3)"].str.contains(r"_", na=False)]
    cl = cl[cl["aaFeatureSequences(CDR3)"].str.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", na=False)]
    cl = cl[cl["readCount"] >= 2]
    print(f"  {len(cl):,} high-quality clones after filtering ({n_before - len(cl):,} removed)")

    cl["V_gene"] = cl["allVHitsWithScore"].apply(clean_gene)
    cl["J_gene"] = cl["allJHitsWithScore"].apply(clean_gene)
    cl["D_gene"] = cl["allDHitsWithScore"].apply(clean_gene)
    cl = cl.dropna(subset=["V_gene", "J_gene"])
    cl["cdr3aa"] = cl["aaFeatureSequences(CDR3)"]
    cl["cdr3nt"] = cl["nFeatureSequences(CDR3)"]
    # Re-normalize freq
    total = cl["readFraction"].sum()
    cl["freq"] = cl["readFraction"] / total
    cl["count"] = cl["readCount"]
    cl["cdr3_len"] = cl["cdr3aa"].str.len()

    # Auto-detect sample ID from clonotype filename
    stem = clonotype_path.stem
    sample_id = stem.split("_gdna_")[0] if "_gdna_" in stem else stem
    # Match TSV row
    tsv_row = tsv[tsv["Sample_ID"].str.startswith(sample_id.split("_")[0])]
    if tsv_row.empty:
        # Try exact match
        tsv_row = tsv[tsv["Sample_ID"].str.contains(stem.split(".")[0][:10])]
    if tsv_row.empty:
        print(f"  WARNING: sample {sample_id} not found in TSV, using mock aggregate data")
        tsv_row = None
    else:
        tsv_row = tsv_row.iloc[0]
        print(f"  Matched TSV row: {tsv_row['Sample_ID']}")

    print(f"Sample ID: {sample_id}")

    # ── 1. sample_info.csv ──────────────────────────────────────────────
    print("Generating CSVs...")
    pd.DataFrame([{
        "SampleID": sample_id,
        "PatientID": sample_id.split("_")[0],
        "Group": "Unknown",
        "Type": "Tissue",
        "SampleType": "gDNA",
    }]).to_csv(output_dir / "sample_info.csv", index=False)

    # ── 2. qc_table.csv (from TSV) ──────────────────────────────────────
    if tsv_row is not None:
        raw_reads = int(tsv_row.get("raw_read", 0))
        clean_reads = int(tsv_row.get("clean_read", 0))
        clonal_reads = int(tsv_row["Clonereads"])
        clean_ratio = round(float(tsv_row["Effective(%)"]) / 100, 4)
        clonal_ratio = round(clonal_reads / clean_reads, 4) if clean_reads > 0 else 0.0
        pd.DataFrame([{
            "Sample": sample_id,
            "rawYield": round(float(tsv_row["Raw_Yield(G)"]), 2),
            "cleanYield": round(float(tsv_row["Clean_Yield(G)"]), 2),
            "rawReads": raw_reads,
            "cleanReads": clean_reads,
            "cleanQ30": round(float(tsv_row["Clean_Q30(%)"]) / 100, 6),
            "cleanRatio": clean_ratio,
            "clonalReads": clonal_reads,
            "clonalRatio": clonal_ratio,
        }]).to_csv(output_dir / "qc_table.csv", index=False)
    else:
        mock_qc(sample_id, output_dir)

    # ── 3. cdr3_length.csv (from clonotype detail) ──────────────────────
    len_freq = cl.groupby("cdr3_len")["freq"].sum()
    len_freq = len_freq / len_freq.sum()
    rows_cdr3 = []
    top_n = 10
    for length in sorted(cl["cdr3_len"].unique()):
        overall = len_freq.get(length, 0)
        rows_cdr3.append({
            "Sample": sample_id, "Length": length,
            "OverallFrequency": round(overall, 6),
            "Clone_ID": "Overall", "CloneFrequency": round(overall, 6),
        })
        df_len = cl[cl["cdr3_len"] == length].nlargest(top_n, "freq")
        for _, r in df_len.iterrows():
            rows_cdr3.append({
                "Sample": sample_id, "Length": length,
                "OverallFrequency": round(overall, 6),
                "Clone_ID": r["cdr3aa"],
                "CloneFrequency": round(r["freq"], 8),
            })
    pd.DataFrame(rows_cdr3).to_csv(output_dir / "cdr3_length.csv", index=False)

    # ── 4. v_gene_freq.csv (from TSV if available, else clonotype) ──────
    if tsv_row is not None:
        v_cols = [c for c in tsv_row.index if isinstance(c, str) and re.match(r"^TRBV\d", c)]
        v_freqs = pd.to_numeric(tsv_row[v_cols], errors="coerce").dropna()
        v_freqs = v_freqs / v_freqs.sum()
        pd.DataFrame({
            "Sample": sample_id, "V_gene": v_freqs.index, "Frequency": v_freqs.values.round(6),
        }).to_csv(output_dir / "v_gene_freq.csv", index=False)
    else:
        v_f = cl.groupby("V_gene")["freq"].sum()
        v_f = v_f / v_f.sum()
        pd.DataFrame({"Sample": sample_id, "V_gene": v_f.index, "Frequency": v_f.values.round(6)}
        ).to_csv(output_dir / "v_gene_freq.csv", index=False)

    # ── 5. j_gene_freq.csv (from TSV) ───────────────────────────────────
    if tsv_row is not None:
        j_cols = [c for c in tsv_row.index if isinstance(c, str) and re.match(r"^TRBJ\d", c)]
        j_freqs = pd.to_numeric(tsv_row[j_cols], errors="coerce").dropna()
        j_freqs = j_freqs / j_freqs.sum()
        pd.DataFrame({
            "Sample": sample_id, "J_gene": j_freqs.index, "Frequency": j_freqs.values.round(6),
        }).to_csv(output_dir / "j_gene_freq.csv", index=False)
    else:
        j_f = cl.groupby("J_gene")["freq"].sum()
        j_f = j_f / j_f.sum()
        pd.DataFrame({"Sample": sample_id, "J_gene": j_f.index, "Frequency": j_f.values.round(6)}
        ).to_csv(output_dir / "j_gene_freq.csv", index=False)

    # ── 6. vj_pairing.csv (from clonotype) ──────────────────────────────
    vj = cl.groupby(["V_gene", "J_gene"])["freq"].sum().reset_index()
    vj["Frequency"] = vj["freq"] / vj["freq"].sum()
    vj = vj[vj["Frequency"] > 0.0001]
    vj["Sample"] = sample_id
    vj[["Sample", "V_gene", "J_gene", "Frequency"]].to_csv(output_dir / "vj_pairing.csv", index=False)

    # ── 7. clone_diversity.csv (from TSV) ───────────────────────────────
    if tsv_row is not None:
        pd.DataFrame([{
            "Sample": sample_id,
            "CloneReads": int(tsv_row["Clonereads"]),
            "SamplingReads": int(tsv_row["SamplingReads"]),
            "CloneCount": int(tsv_row["Clone_AA_count"]),
            "TopFreq": round(float(tsv_row["Top_Freq"]), 6),
            "Clonality": round(float(tsv_row["Clonality_index"]), 6),
            "ShannonIndex": round(float(tsv_row["Shannon_Weiner"]), 6),
            "SimpsonIndex": round(float(tsv_row["Simpson"]), 6),
            "Evenness": round(float(tsv_row["Evenness"]), 6),
        }]).to_csv(output_dir / "clone_diversity.csv", index=False)
    else:
        freqs = cl["freq"].values; freqs = freqs / freqs.sum()
        d = compute_diversity(freqs)
        pd.DataFrame([{
            "Sample": sample_id, "CloneReads": int(cl["readCount"].sum()),
            "SamplingReads": 1_000_000, "CloneCount": len(cl),
            "TopFreq": round(float(freqs.max()), 6), **d,
        }]).to_csv(output_dir / "clone_diversity.csv", index=False)

    # ── 8. clone_composition.csv (from TSV group_freq) ───────────────────
    if tsv_row is not None:
        group_map = {
            "group_freq_2": "Rare (≤0.01%)",
            "group_freq_3": "Small (0.01%-0.1%)",
            "group_freq_4": "Medium (0.1%-1%)",
            "group_freq_5": "Large (1%-10%)",
            "group_freq_6": "Hyperexpanded (>10%)",
            "group_freq_7": "Hyperexpanded (>10%)",
        }
        rows_comp = []
        for col, label in group_map.items():
            v = float(tsv_row.get(col, 0))
            rows_comp.append({"Sample": sample_id, "FrequencyRange": label, "Proportion": round(v, 6)})
        # Merge Hyperexpanded
        merged = {}
        for r in rows_comp:
            k = r["FrequencyRange"]
            merged[k] = merged.get(k, 0) + r["Proportion"]
        final_comp = []
        for r in RANGE_ORDER:
            final_comp.append({"Sample": sample_id, "FrequencyRange": r, "Proportion": round(merged.get(r, 0), 6)})
        pd.DataFrame(final_comp).to_csv(output_dir / "clone_composition.csv", index=False)
    else:
        cl["fr"] = cl["freq"].apply(classify_freq_range)
        comp = cl.groupby("fr")["freq"].count(); comp = comp / comp.sum()
        rows = [{"Sample": sample_id, "FrequencyRange": r, "Proportion": round(float(comp.get(r, 0)), 6)} for r in RANGE_ORDER]
        pd.DataFrame(rows).to_csv(output_dir / "clone_composition.csv", index=False)

    # ── 9. functional_annotation.csv (from clonotype detail + TSV antigen) ──
    # Use clonotype per-clone data; antigen info from TSV is aggregate only
    fa_rows = []
    for _, row in cl.iterrows():
        fa_rows.append({
            "Sample": sample_id, "count": int(row["count"]), "freq": round(row["freq"], 8),
            "cdr3nt": row["cdr3nt"], "cdr3aa": row["cdr3aa"],
            "v": row["V_gene"], "d": row.get("D_gene", ".") or ".",
            "j": row["J_gene"], "VEnd": 0, "DStart": 0, "DEnd": 0, "JStart": 0,
            "id.in.sample": row["cloneId"],
            "match.score": 0.0, "match.weight": 0.0, "gene": "TRB",
            "cdr3": row["cdr3aa"], "species": "HomoSapiens",
            "antigen.epitope": "", "antigen.gene": "",
            "antigen.species": "Unknown",
            "complex.id": "", "v.segm": f"{row['V_gene']}*01",
            "j.segm": f"{row['J_gene']}*01",
            "v.end": 0, "j.start": 0,
            "mhc.a": "", "mhc.b": "", "mhc.class": "",
            "reference.id": "", "vdjdb.score": 0.0,
        })
    pd.DataFrame(fa_rows).to_csv(output_dir / "functional_annotation.csv", index=False)

    # ── 10. cancer_risk.csv (from TSV) ──────────────────────────────────
    if tsv_row is not None:
        pd.DataFrame([{
            "sample_id": sample_id,
            "filtered_records": int(tsv_row.get("filtered_records", len(cl))),
            "KRAS_type": float(tsv_row.get("KRAS_type", 0)),
            "KRAS_score": float(tsv_row.get("KRAS_score", 0)),
            "EGFR_type": float(tsv_row.get("EGFR_type", 0)),
            "EGFR_score": float(tsv_row.get("EGFR_score", 0)),
            "REF_type": float(tsv_row.get("REF_type", 0)),
            "REF_score": float(tsv_row.get("REF_score", 0)),
            "LUNG_CANCER_GDNA_type": float(tsv_row.get("LUNG_CANCER_GDNA_type", 0)),
            "LUNG_CANCER_GDNA_score": float(tsv_row.get("LUNG_CANCER_GDNA_score", 0)),
            "LUNG_CANCER_TISSUE_type": float(tsv_row.get("LUNG_CANCER_TISSUE_type", 0)),
            "LUNG_CANCER_TISSUE_score": float(tsv_row.get("LUNG_CANCER_TISSUE_score", 0)),
        }]).to_csv(output_dir / "cancer_risk.csv", index=False)
    else:
        mock_cancer(sample_id, output_dir, len(cl))

    # ── 10b. antigen_summary.csv (from TSV) ─────────────────────────────
    antigen_rows = []
    seen_names = set()
    if tsv_row is not None:
        # Detailed antigens first (more specific), then summary as fallback
        antigen_cols = [(c, c.replace("_freq", "")) for c in tsv_row.index
                        if isinstance(c, str) and c.endswith("_freq")
                        and not c.startswith("total_freq_")]
        for col, name in antigen_cols:
            f = float(tsv_row.get(col, 0))
            c = int(float(tsv_row.get(f"{name}_count", 0)))
            if f > 0:
                display_name = name.replace("_", " ")
                antigen_rows.append({"Sample": sample_id, "antigen": display_name, "freq": round(f, 8), "count": c})
                seen_names.add(display_name.lower())
        # Summary antigens — only if not already covered by detailed
        summary_antigens = ["CMV", "EBV", "HCV", "HomoSapiens", "InfluenzaA",
                            "SARS-CoV-2", "TriticumAestivum", "YFV"]
        for ag in summary_antigens:
            if ag.lower() in seen_names:
                continue
            f = float(tsv_row.get(f"total_freq_{ag}", 0))
            c = int(float(tsv_row.get(f"count_{ag}", 0)))
            if f > 0:
                antigen_rows.append({"Sample": sample_id, "antigen": ag, "freq": round(f, 8), "count": c})

    pd.DataFrame(antigen_rows).to_csv(output_dir / "antigen_summary.csv", index=False)
    print(f"  Antigen summary: {len(antigen_rows)} antigens")

    # ── 11. tcr_motif.csv (from clonotype) ──────────────────────────────
    top_n = 100
    df_top = cl.nlargest(top_n, "freq")
    motif_rows = []
    for _, row in df_top.iterrows():
        cdr3 = row["cdr3aa"]
        for pos, aa in enumerate(cdr3, 1):
            motif_rows.append({
                "Sample": sample_id,
                "Clone_ID": f"{sample_id}_Clone{row['cloneId']}",
                "CDR3_AA": cdr3, "Position": pos, "AminoAcid": aa,
                "Frequency": round(float(row["freq"]), 6),
            })
    pd.DataFrame(motif_rows).to_csv(output_dir / "tcr_motif.csv", index=False)

    # ── 12. significant_expansion.csv ───────────────────────────────────
    pd.DataFrame(columns=[
        "name", "sample_Ref", "sample_Amp", "fish_pvalue", "p.adj",
        "FC", "count_Ref", "freq_Ref", "count_Amp", "freq_Amp",
    ]).to_csv(output_dir / "significant_expansion.csv", index=False)

    # ── Summary ─────────────────────────────────────────────────────────
    print(f"\nDone! 12 CSVs in {output_dir.resolve()}:")
    for f in sorted(output_dir.glob("*.csv")):
        sz = f.stat().st_size / 1024
        n = pd.read_csv(f).shape[0]
        print(f"  {f.name:40s} {n:>10,} rows  {sz:>8.1f} KB")
    if tsv_row is not None:
        print(f"\nTSV aggregate data used: cancer risk, antigen, diversity, gene freq, QC, composition")
    print(f"Clonotype detail used: CDR3 sequences, V-J pairing, motif, sunburst")


def main():
    p = argparse.ArgumentParser(description="Convert TSV + clonotype to TCR report CSVs")
    p.add_argument("tsv", help="Results_vdj_tcrmatch.tsv")
    p.add_argument("clonotype", help="Clonotype .txt file")
    p.add_argument("-o", "--output-dir", required=True)
    args = p.parse_args()
    convert(args.tsv, args.clonotype, args.output_dir)


if __name__ == "__main__":
    main()
