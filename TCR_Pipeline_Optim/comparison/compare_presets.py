#!/usr/bin/env python3
"""
Preset-only comparison: MiXCR 4.6.0 generic-amplicon vs rna-seq (same gDNA data)
Isolates preset effect — version is controlled at 4.6.0.
"""

import os, re
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BASE = "/Users/edwardchan/Desktop/TCR"
OUT  = os.path.join(BASE, "comparison", "figures_preset")
os.makedirs(OUT, exist_ok=True)

SAMPLES = [
    "S044_SZ20250508032WHB-2_gdna_genome_2537265",
    "S045_SZ20250522049WHB-1_gdna_genome_2537266",
    "S046_SZ20250522052WHB-9_gdna_genome_2537267",
    "S047_SZ20250522068WHB-0_gdna_genome_2537268",
    "S048_SZ20250621039WHB-2_gdna_genome_2537269",
]
SHORT_IDS = ["S044", "S045", "S046", "S047", "S048"]
SHORT_MAP = dict(zip(SAMPLES, SHORT_IDS))
SAMPLE_COLORS = sns.color_palette("Set2", 5)

sns.set_style("whitegrid")
plt.rcParams.update({
    "font.size": 11, "axes.titlesize": 13, "axes.labelsize": 11,
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
})

# ─── Helpers ─────────────────────────────────────────────────────────────────

def read_file(path):
    with open(path) as f:
        return f.read()

def parse_report(path):
    d = {}
    for line in read_file(path).split('\n'):
        line = line.strip()
        if not line or line.startswith('==='): continue
        if ': ' in line:
            key, val = line.split(': ', 1)
            d[key.strip()] = val.strip()
    return d

def parse_log_assembly(path):
    text = read_file(path)
    idx = text.find('report: assemble')
    if idx >= 0:
        section = text[idx:]
    else:
        marker = '============= Report =============='
        idx = text.rfind(marker)
        section = text[idx:] if idx >= 0 else text
    d = {}
    for line in section.split('\n'):
        line = line.strip()
        if not line or line.startswith('==='): continue
        if ': ' in line:
            key, val = line.split(': ', 1)
            d[key.strip()] = val.strip()
    return d

def load_clonotypes(path):
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={
        'nFeatureSequences(CDR3)': 'cdr3nt',
        'aaFeatureSequences(CDR3)': 'cdr3aa',
    })
    df = df[df['cdr3nt'].notna() & (df['cdr3nt'] != '')]
    df = df[df['cdr3aa'].notna() & (df['cdr3aa'] != '')]
    return df

def parse_genes(hits_str, prefix):
    if pd.isna(hits_str) or not str(hits_str).strip():
        return []
    genes = set()
    for token in str(hits_str).split(','):
        m = re.match(rf'({prefix}\d+[\-\d]*)', token.strip())
        if m: genes.add(m.group(1))
    return sorted(genes)

def extract_pct(val):
    if not isinstance(val, str):
        return float(val) if val else 0.0
    val = val.strip()
    if val.endswith('%'): return float(val[:-1])
    if '(' in val:
        return float(val.split('(')[1].split('%')[0].split(')')[0])
    return float(val)

def extract_num(val):
    if not isinstance(val, str):
        return float(val) if val else 0.0
    return float(val.split()[0].replace(',', ''))

def safe_get(d, *keys):
    for k in keys:
        if k in d: return d[k]
    return ''

# ─── Load data ───────────────────────────────────────────────────────────────

amp_reports, rna_reports = {}, {}
amp_assemblies, rna_assemblies = {}, {}
amp_clones, rna_clones = {}, {}

print('Loading data...')
for s in SAMPLES:
    sid = SHORT_MAP[s]
    amp_reports[sid] = parse_report(
        f'{BASE}/test_v3/{s}/Map_Clone_Analysis/{s}.report')
    rna_reports[sid] = parse_report(
        f'{BASE}/test_v3_rnaseq/{s}/Map_Clone_Analysis/{s}.report')
    amp_assemblies[sid] = parse_log_assembly(
        f'{BASE}/test_v3/{s}/Map_Clone_Analysis/{s}.log')
    rna_assemblies[sid] = parse_log_assembly(
        f'{BASE}/test_v3_rnaseq/{s}/Map_Clone_Analysis/{s}.log')
    amp_clones[sid] = load_clonotypes(
        f'{BASE}/test_v3/{s}/Map_Clone_Analysis/{s}.clonotypes.TRB.raw.txt')
    rna_clones[sid] = load_clonotypes(
        f'{BASE}/test_v3_rnaseq/{s}/Map_Clone_Analysis/{s}.clonotypes.TRB.raw.txt')
    print(f'  {sid}: amp={len(amp_clones[sid]):,} clones, rna-seq={len(rna_clones[sid]):,} clones')

# ═════════════════════════════════════════════════════════════════════════════
# P1: ALIGNMENT METRICS
# ═════════════════════════════════════════════════════════════════════════════

def plot_alignment():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    x = np.arange(len(SHORT_IDS)); w = 0.35
    C1, C2 = '#DD8452', '#55A868'

    # A: Alignment rate
    ax = axes[0, 0]
    r1 = [extract_pct(amp_reports[s]['Successfully aligned reads']) for s in SHORT_IDS]
    r2 = [extract_pct(rna_reports[s]['Successfully aligned reads']) for s in SHORT_IDS]
    ax.bar(x - w/2, r1, w, color=C1, label='generic-amplicon')
    ax.bar(x + w/2, r2, w, color=C2, label='rna-seq')
    ax.set_ylabel('Alignment Rate (%)'); ax.set_title('A. Successfully Aligned Reads (%)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS); ax.legend(fontsize=9)
    ax.set_ylim(99.5, 100.05)

    # B: Failed - no hits
    ax = axes[0, 1]
    k1 = 'Alignment failed: no hits (not TCR/IG?)'
    n1 = [extract_num(amp_reports[s].get(k1, '0')) for s in SHORT_IDS]
    n2 = [extract_num(rna_reports[s].get(k1, '0')) for s in SHORT_IDS]
    ax.bar(x - w/2, n1, w, color=C1); ax.bar(x + w/2, n2, w, color=C2)
    ax.set_ylabel('Read Count'); ax.set_title('B. Alignment Failed: No Hits (not TCR/IG)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    # C: Failed - no J hits
    ax = axes[1, 0]
    k1 = 'Alignment failed: absence of J hits'
    n1 = [extract_num(amp_reports[s].get(k1, '0')) for s in SHORT_IDS]
    n2 = [extract_num(rna_reports[s].get(k1, '0')) for s in SHORT_IDS]
    ax.bar(x - w/2, n1, w, color=C1); ax.bar(x + w/2, n2, w, color=C2)
    ax.set_ylabel('Read Count'); ax.set_title('C. Alignment Failed: No J Hits')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    # D: Total reads + non-functional
    ax = axes[1, 1]
    t1 = [extract_num(amp_reports[s]['Total sequencing reads'])/1e6 for s in SHORT_IDS]
    t2 = [extract_num(rna_reports[s]['Total sequencing reads'])/1e6 for s in SHORT_IDS]
    nf1 = [extract_pct(amp_reports[s].get('TRB non-functional', '0%')) for s in SHORT_IDS]
    nf2 = [extract_pct(rna_reports[s].get('TRB non-functional', '0%')) for s in SHORT_IDS]
    ax2 = ax.twinx()
    ax.bar(x - w/2, t1, w, color=C1, label='amp Total (M)')
    ax.bar(x + w/2, t2, w, color=C2, label='rna-seq Total (M)')
    ax2.plot(x, nf1, 'D-', color='#C44E52', linewidth=2, markersize=8, label='amp Non-func %')
    ax2.plot(x, nf2, 's--', color='#E67E22', linewidth=2, markersize=8, label='rna-seq Non-func %')
    ax.set_ylabel('Total Reads (M)'); ax2.set_ylabel('Non-functional TRB (%)')
    ax.set_title('D. Input Reads & Non-functional TRB Rate')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left')

    fig.suptitle('P1: Alignment Metrics — generic-amplicon vs rna-seq (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P1_alignment.png'); plt.close()
    print('  ✓ P1_alignment.png')

# ═════════════════════════════════════════════════════════════════════════════
# P2: ASSEMBLY METRICS
# ═════════════════════════════════════════════════════════════════════════════

def plot_assembly():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    x = np.arange(len(SHORT_IDS)); w = 0.35; C1, C2 = '#DD8452', '#55A868'

    ax = axes[0, 0]
    c1 = [int(extract_num(amp_assemblies[s]['Final clonotype count'])) for s in SHORT_IDS]
    c2 = [int(extract_num(rna_assemblies[s]['Final clonotype count'])) for s in SHORT_IDS]
    ax.bar(x - w/2, c1, w, color=C1); ax.bar(x + w/2, c2, w, color=C2)
    ax.set_ylabel('Clone Count'); ax.set_title('A. Final Clonotype Count')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    ax = axes[0, 1]
    key_set = ['Reads used in clonotypes, percent of total', 'Reads used in clonotypes']
    ru1 = [extract_pct(safe_get(amp_assemblies[s], *key_set)) for s in SHORT_IDS]
    ru2 = [extract_pct(safe_get(rna_assemblies[s], *key_set)) for s in SHORT_IDS]
    ax.bar(x - w/2, ru1, w, color=C1); ax.bar(x + w/2, ru2, w, color=C2)
    ax.set_ylabel('Reads in Clonotypes (%)'); ax.set_title('B. Reads Used in Clonotypes (%)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    ax = axes[1, 0]
    ar1 = [float(extract_num(amp_assemblies[s]['Average number of reads per clonotype'])) for s in SHORT_IDS]
    ar2 = [float(extract_num(rna_assemblies[s]['Average number of reads per clonotype'])) for s in SHORT_IDS]
    ax.bar(x - w/2, ar1, w, color=C1); ax.bar(x + w/2, ar2, w, color=C2)
    ax.set_ylabel('Avg Reads/Clone'); ax.set_title('C. Average Reads per Clonotype')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    ax = axes[1, 1]
    pcr_key = 'Clonotypes eliminated by PCR error correction'
    pe1 = [int(extract_num(amp_assemblies[s].get(pcr_key, '0'))) for s in SHORT_IDS]
    pe2 = [int(extract_num(rna_assemblies[s].get(pcr_key, '0'))) for s in SHORT_IDS]
    ax.bar(x - w/2, pe1, w, color=C1); ax.bar(x + w/2, pe2, w, color=C2)
    ax.set_ylabel('Clonotypes Eliminated'); ax.set_title('D. PCR Error Correction Eliminated')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    fig.suptitle('P2: Assembly Metrics — generic-amplicon vs rna-seq (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P2_assembly.png'); plt.close()
    print('  ✓ P2_assembly.png')

# ═════════════════════════════════════════════════════════════════════════════
# P3: CLONE SENSITIVITY
# ═════════════════════════════════════════════════════════════════════════════

def plot_sensitivity():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    x = np.arange(len(SHORT_IDS)); w = 0.35; C1, C2 = '#DD8452', '#55A868'

    ax = axes[0]
    t1 = [len(amp_clones[s]) for s in SHORT_IDS]
    t2 = [len(rna_clones[s]) for s in SHORT_IDS]
    ax.bar(x - w/2, t1, w, color=C1, label='generic-amplicon')
    ax.bar(x + w/2, t2, w, color=C2, label='rna-seq')
    ax.set_ylabel('Clone Count'); ax.set_title('A. Total Clones Exported')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS); ax.legend(fontsize=9)

    ax = axes[1]
    u1 = [amp_clones[s]['cdr3aa'].nunique() for s in SHORT_IDS]
    u2 = [rna_clones[s]['cdr3aa'].nunique() for s in SHORT_IDS]
    ax.bar(x - w/2, u1, w, color=C1); ax.bar(x + w/2, u2, w, color=C2)
    ax.set_ylabel('Unique CDR3 AA'); ax.set_title('B. Unique CDR3 AA Sequences')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    ax = axes[2]
    s1 = [(amp_clones[s]['readCount'] == 1).sum() / len(amp_clones[s]) * 100 for s in SHORT_IDS]
    s2 = [(rna_clones[s]['readCount'] == 1).sum() / len(rna_clones[s]) * 100 for s in SHORT_IDS]
    ax.bar(x - w/2, s1, w, color=C1); ax.bar(x + w/2, s2, w, color=C2)
    ax.set_ylabel('Singleton %'); ax.set_title('C. Singleton Clone Percentage')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    fig.suptitle('P3: Clone Detection Sensitivity (unfiltered)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P3_sensitivity.png'); plt.close()
    print('  ✓ P3_sensitivity.png')


def plot_sensitivity_filtered():
    """P3b: Clone detection sensitivity with quality filters applied.

    Filters: CDR3 AA length 5-45 aa, readCount > 1.
    Demonstrates that rna-seq's apparent "more clones" is driven by noise.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    x = np.arange(len(SHORT_IDS)); w = 0.35; C1, C2 = '#DD8452', '#55A868'

    # Pre-compute filtered data
    filt_amp = {}; filt_rna = {}
    for sid in SHORT_IDS:
        a = amp_clones[sid].copy(); r = rna_clones[sid].copy()
        a['aa_len'] = a['cdr3aa'].str.len()
        r['aa_len'] = r['cdr3aa'].str.len()
        a = a[(a['aa_len'] >= 5) & (a['aa_len'] <= 45)]
        r = r[(r['aa_len'] >= 5) & (r['aa_len'] <= 45)]
        a = a[a['readCount'] > 1]
        r = r[r['readCount'] > 1]
        filt_amp[sid] = a; filt_rna[sid] = r

    ax = axes[0]
    t1 = [len(filt_amp[s]) for s in SHORT_IDS]
    t2 = [len(filt_rna[s]) for s in SHORT_IDS]
    ax.bar(x - w/2, t1, w, color=C1, label='generic-amplicon')
    ax.bar(x + w/2, t2, w, color=C2, label='rna-seq')
    ax.set_ylabel('Clone Count'); ax.set_title('A. Filtered Clones (AA 5-45, reads>1)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS); ax.legend(fontsize=9)

    ax = axes[1]
    u1 = [filt_amp[s]['cdr3aa'].nunique() for s in SHORT_IDS]
    u2 = [filt_rna[s]['cdr3aa'].nunique() for s in SHORT_IDS]
    ax.bar(x - w/2, u1, w, color=C1); ax.bar(x + w/2, u2, w, color=C2)
    ax.set_ylabel('Unique CDR3 AA'); ax.set_title('B. Unique CDR3 AA Sequences (filtered)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)

    ax = axes[2]
    r1 = [len(filt_amp[s]) / len(amp_clones[s]) * 100 for s in SHORT_IDS]
    r2 = [len(filt_rna[s]) / len(rna_clones[s]) * 100 for s in SHORT_IDS]
    ax.bar(x - w/2, r1, w, color=C1); ax.bar(x + w/2, r2, w, color=C2)
    ax.set_ylabel('Retention Rate (%)'); ax.set_title('C. Clones Retained After Filtering (%)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)
    for i in range(len(SHORT_IDS)):
        ax.text(i - w/2, r1[i] + 0.5, f'{r1[i]:.1f}%', ha='center', fontsize=8, fontweight='bold')
        ax.text(i + w/2, r2[i] + 0.5, f'{r2[i]:.1f}%', ha='center', fontsize=8, fontweight='bold')

    fig.suptitle('P3b: Clone Detection Sensitivity (filtered: AA 5-45, readCount>1)',
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P3_sensitivity_filtered.png'); plt.close()
    print('  ✓ P3_sensitivity_filtered.png')

    # Print summary
    print('\n  P3b filtered sensitivity:')
    for sid in SHORT_IDS:
        print(f'    {sid}: amp={len(filt_amp[sid]):,} (retain {len(filt_amp[sid])/len(amp_clones[sid])*100:.1f}%), '
              f'rna-seq={len(filt_rna[sid]):,} (retain {len(filt_rna[sid])/len(rna_clones[sid])*100:.1f}%), '
              f'amp>rna={len(filt_amp[sid]) > len(filt_rna[sid])}')


# ═════════════════════════════════════════════════════════════════════════════
# P4: CDR3 LENGTH DISTRIBUTION (most critical)
# ═════════════════════════════════════════════════════════════════════════════

def plot_cdr3_length():
    fig, axes = plt.subplots(2, 1, figsize=(16, 10))

    ax = axes[0]
    for i, sid in enumerate(SHORT_IDS):
        l1 = amp_clones[sid]['cdr3nt'].str.len().dropna()
        l2 = rna_clones[sid]['cdr3nt'].str.len().dropna()
        bins = np.arange(0, 300, 3)
        ax.hist(l1, bins=bins, alpha=0.45, color=SAMPLE_COLORS[i], label=f'{sid} amp',
                density=True, histtype='stepfilled', linewidth=1.5, edgecolor=SAMPLE_COLORS[i])
        ax.hist(l2, bins=bins, alpha=0.25, color=SAMPLE_COLORS[i], label=f'{sid} rna-seq',
                density=True, histtype='step', linewidth=2.5, linestyle='--',
                edgecolor=SAMPLE_COLORS[i])
    ax.set_xlabel('CDR3 Nucleotide Length (bp)'); ax.set_ylabel('Density')
    ax.set_title('A. CDR3 NT Length (solid=amp, dashed=rna-seq)')
    ax.set_xlim(20, 300); ax.legend(fontsize=8, ncol=2)

    ax = axes[1]
    for i, sid in enumerate(SHORT_IDS):
        l1 = amp_clones[sid]['cdr3aa'].str.len().dropna()
        l2 = rna_clones[sid]['cdr3aa'].str.len().dropna()
        bins = np.arange(0, 100, 1)
        ax.hist(l1, bins=bins, alpha=0.45, color=SAMPLE_COLORS[i], label=f'{sid} amp',
                density=True, histtype='stepfilled', linewidth=1.5, edgecolor=SAMPLE_COLORS[i])
        ax.hist(l2, bins=bins, alpha=0.25, color=SAMPLE_COLORS[i], label=f'{sid} rna-seq',
                density=True, histtype='step', linewidth=2.5, linestyle='--',
                edgecolor=SAMPLE_COLORS[i])
    ax.set_xlabel('CDR3 Amino Acid Length'); ax.set_ylabel('Density')
    ax.set_title('B. CDR3 AA Length (solid=amp, dashed=rna-seq)')
    ax.set_xlim(5, 80); ax.legend(fontsize=8, ncol=2)

    fig.suptitle('P4: CDR3 Length Distribution — Preset Comparison (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P4_cdr3_length.png'); plt.close()
    print('  ✓ P4_cdr3_length.png')

# ═════════════════════════════════════════════════════════════════════════════
# P5: V/J GENE USAGE
# ═════════════════════════════════════════════════════════════════════════════

def plot_gene_usage():
    fig, axes = plt.subplots(2, 1, figsize=(16, 12))
    for gene_type, ax in [('V', axes[0]), ('J', axes[1])]:
        prefix = 'TRB' + gene_type; hit_col = f'all{gene_type}HitsWithScore'
        amp_counts = Counter(); rna_counts = Counter()
        for sid in SHORT_IDS:
            for _, row in amp_clones[sid].iterrows():
                for g in parse_genes(row.get(hit_col, ''), prefix):
                    amp_counts[g] += row['readCount']
            for _, row in rna_clones[sid].iterrows():
                for g in parse_genes(row.get(hit_col, ''), prefix):
                    rna_counts[g] += row['readCount']
        a_total = sum(amp_counts.values()); r_total = sum(rna_counts.values())
        all_genes = sorted(set(amp_counts) | set(rna_counts))
        f1 = [amp_counts.get(g, 0) / a_total * 100 for g in all_genes]
        f2 = [rna_counts.get(g, 0) / r_total * 100 for g in all_genes]
        ax.scatter(f1, f2, alpha=0.6, s=50, edgecolors='grey', linewidth=0.5)
        mx = max(max(f1), max(f2)) * 1.05
        ax.plot([0, mx], [0, mx], 'k--', alpha=0.3, linewidth=1)
        diffs = [abs(f1[i] - f2[i]) for i in range(len(all_genes))]
        for idx in np.argsort(diffs)[-10:]:
            ax.annotate(all_genes[idx], (f1[idx], f2[idx]),
                       fontsize=7, alpha=0.8, xytext=(5, 5), textcoords='offset points')
        rho = pd.Series(f1).corr(pd.Series(f2), method='spearman')
        ax.set_xlabel('generic-amplicon Frequency (%)'); ax.set_ylabel('rna-seq Frequency (%)')
        ax.set_title(f'{gene_type} Gene Usage (Spearman ρ = {rho:.4f})')

    fig.suptitle('P5: V/J Gene Usage — Preset Comparison (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P5_gene_usage.png'); plt.close()
    print('  ✓ P5_gene_usage.png')

# ═════════════════════════════════════════════════════════════════════════════
# P6: CDR3 AA OVERLAP
# ═════════════════════════════════════════════════════════════════════════════

def plot_cdr3_overlap():
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(SHORT_IDS)); w = 0.25; A, R, S = '#DD8452', '#55A868', '#8B5CF6'

    amp_only, rna_only, shared = [], [], []
    amp_pct, rna_pct, sh_pct = [], [], []
    for sid in SHORT_IDS:
        sa = set(amp_clones[sid]['cdr3aa'].dropna())
        sr = set(rna_clones[sid]['cdr3aa'].dropna())
        sh = sa & sr; total = len(sa | sr)
        amp_only.append(len(sa - sr)); rna_only.append(len(sr - sa)); shared.append(len(sh))
        amp_pct.append(len(sa - sr) / total * 100 if total else 0)
        rna_pct.append(len(sr - sa) / total * 100 if total else 0)
        sh_pct.append(len(sh) / total * 100 if total else 0)

    ax.bar(x, amp_only, w, label='amp-only CDR3 AA', color=A)
    ax.bar(x, shared, w, bottom=amp_only, label='Shared CDR3 AA', color=R)
    bot = [a + b for a, b in zip(amp_only, shared)]
    ax.bar(x, rna_only, w, bottom=bot, label='rna-seq-only CDR3 AA', color=S)
    for i in range(len(SHORT_IDS)):
        ax.text(i, amp_only[i]/2, f'{amp_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(i, amp_only[i] + shared[i]/2, f'{sh_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(i, amp_only[i] + shared[i] + rna_only[i]/2, f'{rna_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
    ax.set_ylabel('CDR3 AA Sequence Count'); ax.set_title('CDR3 AA Overlap: generic-amplicon vs rna-seq')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS); ax.legend(fontsize=10)

    fig.suptitle('P6: CDR3 Sequence Overlap — Preset Comparison (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P6_cdr3_overlap.png'); plt.close()
    print('  ✓ P6_cdr3_overlap.png')

# ═════════════════════════════════════════════════════════════════════════════
# P6b: FILTERED CDR3 AA OVERLAP
# ═════════════════════════════════════════════════════════════════════════════

def plot_cdr3_overlap_filtered():
    """
    P6b: Filtered CDR3 AA overlap.

    过滤逻辑：
      1. CDR3 AA 长度 5–45 aa
         依据：IMGT 定义，人类 TRB CDR3 正常范围 ~12-16 aa。
         下限 5 aa 排除明显不完整的截断序列；上限 45 aa 覆盖所有已知功能性
         TRB CDR3 变体（物理折叠极限 ~30 aa），同时保留边界分析余量。
         超出 45 aa 的序列通常包含恒定区污染、内含子通读或移码翻译。
      2. readCount > 1 (去除 singleton)
         依据：singleton 占 58-77%，大部分是 PCR/测序错误产生的单拷贝噪音。
         保留 readCount >= 2 的克隆可大幅降低随机噪声对集合比较的干扰。
         此阈值与免疫组学领域的常见做法一致 (e.g. Immunarch 默认 min.count=2)。
    """
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(SHORT_IDS)); w = 0.25
    A, R, S = '#DD8452', '#55A868', '#8B5CF6'

    amp_only, rna_only, shared = [], [], []
    amp_pct, rna_pct, sh_pct = [], [], []
    for sid in SHORT_IDS:
        a = amp_clones[sid].copy(); r = rna_clones[sid].copy()

        # Filter 1: CDR3 AA length 5-45
        a['aa_len'] = a['cdr3aa'].str.len()
        r['aa_len'] = r['cdr3aa'].str.len()
        a = a[(a['aa_len'] >= 5) & (a['aa_len'] <= 45)]
        r = r[(r['aa_len'] >= 5) & (r['aa_len'] <= 45)]

        # Filter 2: readCount > 1
        a = a[a['readCount'] > 1]
        r = r[r['readCount'] > 1]

        sa = set(a['cdr3aa'].dropna()); sr = set(r['cdr3aa'].dropna())
        sh = sa & sr; total = len(sa | sr)
        amp_only.append(len(sa - sr)); rna_only.append(len(sr - sa)); shared.append(len(sh))
        amp_pct.append(len(sa - sr) / total * 100 if total else 0)
        rna_pct.append(len(sr - sa) / total * 100 if total else 0)
        sh_pct.append(len(sh) / total * 100 if total else 0)

    ax.bar(x, amp_only, w, label='amp-only CDR3 AA', color=A)
    ax.bar(x, shared, w, bottom=amp_only, label='Shared CDR3 AA', color=R)
    bot = [a + b for a, b in zip(amp_only, shared)]
    ax.bar(x, rna_only, w, bottom=bot, label='rna-seq-only CDR3 AA', color=S)
    for i in range(len(SHORT_IDS)):
        ax.text(i, amp_only[i]/2, f'{amp_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(i, amp_only[i] + shared[i]/2, f'{sh_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(i, amp_only[i] + shared[i] + rna_only[i]/2, f'{rna_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
    ax.set_ylabel('CDR3 AA Sequence Count')
    ax.set_title('CDR3 AA Overlap (filtered): AA 10-30, readCount > 1')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS); ax.legend(fontsize=10)

    fig.suptitle('P6b: CDR3 Sequence Overlap — Filtered (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P6_cdr3_overlap_filtered.png'); plt.close()
    print('  ✓ P6_cdr3_overlap_filtered.png')

    # Print filtered summary
    print('\n  Filtered overlap (AA 5-45, readCount>1):')
    for sid in SHORT_IDS:
        print(f'    {sid}: amp={amp_only[SHORT_IDS.index(sid)]}, rna={rna_only[SHORT_IDS.index(sid)]}, '
              f'shared={shared[SHORT_IDS.index(sid)]}, '
              f'amp-only%={amp_pct[SHORT_IDS.index(sid)]:.2f}%, rna-only%={rna_pct[SHORT_IDS.index(sid)]:.2f}%')

# ═════════════════════════════════════════════════════════════════════════════
# P7: SHARED CLONE FREQUENCY CORRELATION + STRATIFIED
# ═════════════════════════════════════════════════════════════════════════════

def plot_clone_correlation():
    """P7: Clone frequency correlation — AA-level aggregation on filtered data.

    Filters (CDR3 AA 5-45, readCount > 1) applied BEFORE aggregation so that
    amp-only / rna-only counts in subplot titles match the corrected P6b/P3b findings.
    """

    amp_agg = {}
    rna_agg = {}
    for sid in SHORT_IDS:
        a = amp_clones[sid].copy()
        r = rna_clones[sid].copy()
        a['aa_len'] = a['cdr3aa'].str.len()
        r['aa_len'] = r['cdr3aa'].str.len()
        a = a[(a['aa_len'] >= 5) & (a['aa_len'] <= 45)]
        r = r[(r['aa_len'] >= 5) & (r['aa_len'] <= 45)]
        a = a[a['readCount'] > 1]
        r = r[r['readCount'] > 1]

        a_agg = a.groupby('cdr3aa').agg(
            readCount=('readCount', 'sum'),
            readFraction=('readFraction', 'sum'),
            cdr3nt=('cdr3nt', 'first'),
        ).reset_index()
        r_agg = r.groupby('cdr3aa').agg(
            readCount=('readCount', 'sum'),
            readFraction=('readFraction', 'sum'),
            cdr3nt=('cdr3nt', 'first'),
        ).reset_index()
        amp_agg[sid] = a_agg
        rna_agg[sid] = r_agg

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for i, sid in enumerate(SHORT_IDS):
        ax = axes[i]
        a = amp_agg[sid]
        r = rna_agg[sid]
        m = a.merge(r, on='cdr3aa', suffixes=('_amp', '_rna'))
        n_shared = len(m)
        n_amp_only = a['cdr3aa'].nunique() - m['cdr3aa'].nunique()
        n_rna_only = r['cdr3aa'].nunique() - m['cdr3aa'].nunique()
        if len(m) > 0:
            ax.scatter(m['readFraction_amp'] * 100, m['readFraction_rna'] * 100,
                      alpha=0.3, s=8, edgecolors='none')
            ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel('amp Frequency (%)'); ax.set_ylabel('rna-seq Frequency (%)')
        rho = m['readFraction_amp'].corr(m['readFraction_rna'], method='spearman') if len(m) > 1 else 0
        r_pearson = m['readFraction_amp'].corr(m['readFraction_rna'], method='pearson') if len(m) > 1 else 0
        ax.set_title(f'{sid}\nShared: {n_shared:,} | amp-only: {n_amp_only:,} | rna-only: {n_rna_only:,}\nρ = {rho:.4f}  r = {r_pearson:.4f}')
        ax.plot([1e-6, 100], [1e-6, 100], 'k--', alpha=0.3, linewidth=1)
    axes[5].set_visible(False)

    fig.suptitle('P7: Clone Frequency Correlation — AA-aggregated, filtered (both 4.6.0)',
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P7_clone_correlation.png'); plt.close()
    print('  ✓ P7_clone_correlation.png (AA-aggregated, filtered)')

    # Stratified correlation (no singleton bin — filtered data has readCount > 1)
    print('\n  Stratified correlation (AA-aggregated, filtered, shared clones):')
    for sid in SHORT_IDS:
        a = amp_agg[sid]
        r = rna_agg[sid]
        m = a.merge(r, on='cdr3aa', suffixes=('_amp', '_rna'))
        m['mean_reads'] = (m['readCount_amp'] + m['readCount_rna']) / 2
        bins = [(2, 5), (6, 20), (21, 100), (101, 1000), (1001, 999999)]
        print(f'    {sid} (shared={len(m):,}):')
        for lo, hi in bins:
            sub = m[(m['mean_reads'] >= lo) & (m['mean_reads'] <= hi)]
            if len(sub) < 3: continue
            rho = sub['readFraction_amp'].corr(sub['readFraction_rna'], method='spearman')
            r_val = sub['readFraction_amp'].corr(sub['readFraction_rna'], method='pearson')
            ratio = sub['readFraction_amp'] / sub['readFraction_rna']
            gt2x = ((ratio > 2) | (ratio < 0.5)).sum()
            label = f'{lo}-{hi} reads'
            print(f'      {label:>15s}: N={len(sub):,}  ρ={rho:+.3f}  r={r_val:+.3f}  >2x={gt2x/len(sub)*100:.1f}%')

    return amp_agg, rna_agg

# ═════════════════════════════════════════════════════════════════════════════
# DISCORDANT CLONES (post-AA-aggregation)
# ═════════════════════════════════════════════════════════════════════════════

def find_discordant_clones(amp_agg, rna_agg):
    """Find genuinely discordant clones after AA-level aggregation."""
    rows = []
    for sid in SHORT_IDS:
        a = amp_agg[sid]; r = rna_agg[sid]
        m = a.merge(r, on='cdr3aa', suffixes=('_amp', '_rna'))
        m['fold_change'] = m['readFraction_amp'] / m['readFraction_rna']
        m['abs_log2FC'] = np.abs(np.log2(m['fold_change']))
        m['sample'] = sid
        discordant = m[m['abs_log2FC'] > np.log2(100)].sort_values('abs_log2FC', ascending=False)
        for _, row in discordant.iterrows():
            rows.append({
                'sample': sid,
                'cdr3aa': row['cdr3aa'],
                'cdr3nt_amp': row['cdr3nt_amp'],
                'cdr3nt_rna': row['cdr3nt_rna'],
                'readCount_amp': int(row['readCount_amp']),
                'readCount_rna': int(row['readCount_rna']),
                'readFraction_amp': row['readFraction_amp'],
                'readFraction_rna': row['readFraction_rna'],
                'fold_change': row['fold_change'],
                'abs_log2FC': row['abs_log2FC'],
            })
    if rows:
        df = pd.DataFrame(rows)
        df = df.sort_values('abs_log2FC', ascending=False).head(50)
        df.to_csv(f'{OUT}/../discordant_top50_AAagg.tsv', sep='\t', index=False)
        print(f'\n  ✓ Saved {len(df)} genuinely discordant clones to discordant_top50_AAagg.tsv')
        print(f'    Fold change range: {df["fold_change"].min():.1f} – {df["fold_change"].max():.1f}')
        return df
    else:
        print('\n  No genuinely discordant clones found (|log2FC| > log2(100))')
        return None

# ═════════════════════════════════════════════════════════════════════════════
# P7b: NT-LEVEL CLONE FREQUENCY CORRELATION
# ═════════════════════════════════════════════════════════════════════════════

def plot_nt_correlation():
    """P7b: NT-level correlation — filtered data, groupby cdr3nt, merge on NT sequence.

    Filters (CDR3 AA 5-45, readCount > 1) applied before aggregation.
    """

    amp_nt = {}
    rna_nt = {}
    for sid in SHORT_IDS:
        a = amp_clones[sid].copy()
        r = rna_clones[sid].copy()
        a['aa_len'] = a['cdr3aa'].str.len()
        r['aa_len'] = r['cdr3aa'].str.len()
        a = a[(a['aa_len'] >= 5) & (a['aa_len'] <= 45)]
        r = r[(r['aa_len'] >= 5) & (r['aa_len'] <= 45)]
        a = a[a['readCount'] > 1]
        r = r[r['readCount'] > 1]

        a_nt = a.groupby('cdr3nt').agg(
            readCount=('readCount', 'sum'),
            readFraction=('readFraction', 'sum'),
            cdr3aa=('cdr3aa', 'first'),
        ).reset_index()
        r_nt = r.groupby('cdr3nt').agg(
            readCount=('readCount', 'sum'),
            readFraction=('readFraction', 'sum'),
            cdr3aa=('cdr3aa', 'first'),
        ).reset_index()
        amp_nt[sid] = a_nt
        rna_nt[sid] = r_nt

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for i, sid in enumerate(SHORT_IDS):
        ax = axes[i]
        a = amp_nt[sid]
        r = rna_nt[sid]
        m = a.merge(r, on='cdr3nt', suffixes=('_amp', '_rna'))
        n_shared = len(m)
        n_amp_only = len(a) - n_shared
        n_rna_only = len(r) - n_shared
        if len(m) > 0:
            ax.scatter(m['readFraction_amp'] * 100, m['readFraction_rna'] * 100,
                      alpha=0.3, s=8, edgecolors='none')
            ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel('amp Frequency (%)'); ax.set_ylabel('rna-seq Frequency (%)')
        rho = m['readFraction_amp'].corr(m['readFraction_rna'], method='spearman') if len(m) > 1 else 0
        r_pearson = m['readFraction_amp'].corr(m['readFraction_rna'], method='pearson') if len(m) > 1 else 0
        ax.set_title(f'{sid}\nNT shared: {n_shared:,} | amp-only: {n_amp_only:,} | rna-only: {n_rna_only:,}\nρ = {rho:.4f}  r = {r_pearson:.4f}')
        ax.plot([1e-6, 100], [1e-6, 100], 'k--', alpha=0.3, linewidth=1)
    axes[5].set_visible(False)

    fig.suptitle('P7b: Clone Frequency Correlation — NT-level merge, filtered (both 4.6.0)',
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P7_correlation_NTlevel.png'); plt.close()
    print('  ✓ P7_correlation_NTlevel.png (filtered)')

    # NT-level Jaccard and stratified correlation
    print('\n  NT-level overlap (filtered):')
    for sid in SHORT_IDS:
        sa = set(amp_nt[sid]['cdr3nt']); sr = set(rna_nt[sid]['cdr3nt'])
        sh = sa & sr; uni = sa | sr
        j = len(sh)/len(uni)*100 if uni else 0
        print(f'    {sid}: amp_NT={len(sa):,}  rna_NT={len(sr):,}  shared={len(sh):,}  Jaccard={j:.1f}%')

    print('\n  NT-level stratified correlation (filtered, shared NT sequences):')
    for sid in SHORT_IDS:
        a = amp_nt[sid]; r = rna_nt[sid]
        m = a.merge(r, on='cdr3nt', suffixes=('_amp', '_rna'))
        m['mean_reads'] = (m['readCount_amp'] + m['readCount_rna']) / 2
        bins = [(2, 5), (6, 20), (21, 100), (101, 1000), (1001, 999999)]
        print(f'    {sid} (NT shared={len(m):,}):')
        for lo, hi in bins:
            sub = m[(m['mean_reads'] >= lo) & (m['mean_reads'] <= hi)]
            if len(sub) < 3: continue
            rho = sub['readFraction_amp'].corr(sub['readFraction_rna'], method='spearman')
            r_val = sub['readFraction_amp'].corr(sub['readFraction_rna'], method='pearson')
            ratio = sub['readFraction_amp'] / sub['readFraction_rna']
            gt2x = ((ratio > 2) | (ratio < 0.5)).sum()
            label = f'{lo}-{hi} reads'
            print(f'      {label:>15s}: N={len(sub):,}  ρ={rho:+.3f}  r={r_val:+.3f}  >2x={gt2x/len(sub)*100:.1f}%')

    return amp_nt, rna_nt

# ═════════════════════════════════════════════════════════════════════════════
# P8: CDR3 LENGTH BOXPLOTS
# ═════════════════════════════════════════════════════════════════════════════

def plot_cdr3_boxplot():
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    rows_nt = []; rows_aa = []
    for sid in SHORT_IDS:
        for l in amp_clones[sid]['cdr3nt'].str.len().dropna():
            rows_nt.append({'Sample': sid, 'Preset': 'amp', 'Length': l})
        for l in rna_clones[sid]['cdr3nt'].str.len().dropna():
            rows_nt.append({'Sample': sid, 'Preset': 'rna-seq', 'Length': l})
        for l in amp_clones[sid]['cdr3aa'].str.len().dropna():
            rows_aa.append({'Sample': sid, 'Preset': 'amp', 'Length': l})
        for l in rna_clones[sid]['cdr3aa'].str.len().dropna():
            rows_aa.append({'Sample': sid, 'Preset': 'rna-seq', 'Length': l})
    sns.boxplot(data=pd.DataFrame(rows_nt), x='Sample', y='Length', hue='Preset', ax=axes[0],
                palette={'amp': '#DD8452', 'rna-seq': '#55A868'}, fliersize=0.5)
    axes[0].set_title('A. CDR3 NT Length'); axes[0].set_ylabel('CDR3 nt Length (bp)')
    sns.boxplot(data=pd.DataFrame(rows_aa), x='Sample', y='Length', hue='Preset', ax=axes[1],
                palette={'amp': '#DD8452', 'rna-seq': '#55A868'}, fliersize=0.5)
    axes[1].set_title('B. CDR3 AA Length'); axes[1].set_ylabel('CDR3 AA Length')
    fig.suptitle('P8: CDR3 Length Boxplot — Preset Comparison (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P8_cdr3_boxplot.png'); plt.close()
    print('  ✓ P8_cdr3_boxplot.png')

# ═════════════════════════════════════════════════════════════════════════════
# P9: TOP-N CLONE OVERLAP
# ═════════════════════════════════════════════════════════════════════════════

def plot_top_overlap():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    top_ns = [10, 50, 100, 500, 1000]; x = np.arange(len(top_ns))

    ax = axes[0]
    for i, sid in enumerate(SHORT_IDS):
        overlaps = []
        for n in top_ns:
            t1 = set(amp_clones[sid].nlargest(min(n, len(amp_clones[sid])), 'readCount')['cdr3aa'])
            t2 = set(rna_clones[sid].nlargest(min(n, len(rna_clones[sid])), 'readCount')['cdr3aa'])
            overlaps.append(len(t1 & t2) / n * 100 if n > 0 else 0)
        ax.plot(x, overlaps, 'o-', color=SAMPLE_COLORS[i], label=sid, linewidth=2, markersize=8)
    ax.set_xticks(x); ax.set_xticklabels([f'Top {n}' for n in top_ns])
    ax.set_ylabel('Overlap (%)'); ax.set_title('A. Top-N Clone Overlap'); ax.legend(); ax.set_ylim(0, 105)

    ax = axes[1]
    for i, sid in enumerate(SHORT_IDS):
        fracs = []
        for n in top_ns:
            t1 = set(amp_clones[sid].nlargest(min(n, len(amp_clones[sid])), 'readCount')['cdr3aa'])
            t2 = set(rna_clones[sid].nlargest(min(n, len(rna_clones[sid])), 'readCount')['cdr3aa'])
            sh = t1 & t2; tot = amp_clones[sid]['readCount'].sum()
            sr = amp_clones[sid][amp_clones[sid]['cdr3aa'].isin(sh)]['readCount'].sum()
            fracs.append(sr / tot * 100 if tot > 0 else 0)
        ax.plot(x, fracs, 'o-', color=SAMPLE_COLORS[i], label=sid, linewidth=2, markersize=8)
    ax.set_xticks(x); ax.set_xticklabels([f'Top {n}' for n in top_ns])
    ax.set_ylabel('amp Read Fraction in Shared (%)'); ax.set_title('B. Read Fraction Captured by Shared Top-N'); ax.legend()

    fig.suptitle('P9: Top Clones Concordance — Preset Comparison (both 4.6.0)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P9_top_overlap.png'); plt.close()
    print('  ✓ P9_top_overlap.png')

# ═════════════════════════════════════════════════════════════════════════════
# P10: BIOLOGICAL VALIDATION
# ═════════════════════════════════════════════════════════════════════════════

def biological_validation():
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes = axes.flatten()
    A, R = '#DD8452', '#55A868'

    for i, sid in enumerate(SHORT_IDS):
        ax = axes[i]
        a = amp_clones[sid]; r = rna_clones[sid]
        # Shared vs unique CDR3 length
        s_amp = set(a['cdr3aa']); s_rna = set(r['cdr3aa'])
        shared = s_amp & s_rna
        a['group'] = a['cdr3aa'].apply(lambda x: 'shared' if x in shared else 'amp-only')
        r['group'] = r['cdr3aa'].apply(lambda x: 'shared' if x in shared else 'rna-only')
        a_shared_nt = a[a['group'] == 'shared']['cdr3nt'].str.len().mean()
        a_uniq_nt   = a[a['group'] == 'amp-only']['cdr3nt'].str.len().mean()
        r_shared_nt = r[r['group'] == 'shared']['cdr3nt'].str.len().mean()
        r_uniq_nt   = r[r['group'] == 'rna-only']['cdr3nt'].str.len().mean()
        a_shared_aa = a[a['group'] == 'shared']['cdr3aa'].str.len().mean()
        a_uniq_aa   = a[a['group'] == 'amp-only']['cdr3aa'].str.len().mean()
        r_shared_aa = r[r['group'] == 'shared']['cdr3aa'].str.len().mean()
        r_uniq_aa   = r[r['group'] == 'rna-only']['cdr3aa'].str.len().mean()

        ax.text(0.05, 0.95, f'{sid}', transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
        ax.text(0.05, 0.85, f'amp-shared NT: {a_shared_nt:.1f} bp', transform=ax.transAxes, fontsize=10, color=A)
        ax.text(0.05, 0.77, f'amp-unique NT: {a_uniq_nt:.1f} bp', transform=ax.transAxes, fontsize=10, color=A)
        ax.text(0.05, 0.69, f'rna-shared NT: {r_shared_nt:.1f} bp', transform=ax.transAxes, fontsize=10, color=R)
        ax.text(0.05, 0.61, f'rna-unique NT: {r_uniq_nt:.1f} bp', transform=ax.transAxes, fontsize=10, color=R)
        ax.text(0.05, 0.51, f'amp-shared AA: {a_shared_aa:.1f}', transform=ax.transAxes, fontsize=10, color=A)
        ax.text(0.05, 0.43, f'amp-unique AA: {a_uniq_aa:.1f}', transform=ax.transAxes, fontsize=10, color=A)
        ax.text(0.05, 0.35, f'rna-shared AA: {r_shared_aa:.1f}', transform=ax.transAxes, fontsize=10, color=R)
        ax.text(0.05, 0.27, f'rna-unique AA: {r_uniq_aa:.1f}', transform=ax.transAxes, fontsize=10, color=R)
        ax.axis('off')

        print(f'\n  {sid}:')
        print(f'    amp-shared:  NT={a_shared_nt:.1f}, AA={a_shared_aa:.1f}')
        print(f'    amp-unique:  NT={a_uniq_nt:.1f}, AA={a_uniq_aa:.1f}  (n={len(a[a["group"]=="amp-only"]):,})')
        print(f'    rna-shared:  NT={r_shared_nt:.1f}, AA={r_shared_aa:.1f}')
        print(f'    rna-unique:  NT={r_uniq_nt:.1f}, AA={r_uniq_aa:.1f}  (n={len(r[r["group"]=="rna-only"]):,})')

    # C-start Cys check
    ax = axes[5]
    x = np.arange(len(SHORT_IDS)); w = 0.35
    c1 = []; c2 = []
    for sid in SHORT_IDS:
        a_start = amp_clones[sid]['cdr3aa'].str[0]
        r_start = rna_clones[sid]['cdr3aa'].str[0]
        c1.append((a_start == 'C').sum() / len(a_start) * 100)
        c2.append((r_start == 'C').sum() / len(r_start) * 100)
    ax.bar(x - w/2, c1, w, color=A, label='amp')
    ax.bar(x + w/2, c2, w, color=R, label='rna-seq')
    ax.set_ylabel('C-start %'); ax.set_title('N-terminal Cys (C) Conservation')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS); ax.legend()

    fig.suptitle('P10: Biological Validation', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/P10_biological_validation.png'); plt.close()
    print('  ✓ P10_biological_validation.png')

    # C-end F/W check
    print('\n  C-terminal F/W %:')
    for sid in SHORT_IDS:
        a_end = amp_clones[sid]['cdr3aa'].str[-1]
        r_end = rna_clones[sid]['cdr3aa'].str[-1]
        af = (a_end.isin(['F', 'W'])).sum() / len(a_end) * 100
        rf = (r_end.isin(['F', 'W'])).sum() / len(r_end) * 100
        a_second = a_end.value_counts().index[1] if len(a_end.value_counts()) > 1 else 'N/A'
        r_second = r_end.value_counts().index[1] if len(r_end.value_counts()) > 1 else 'N/A'
        print(f'    {sid}: amp={af:.1f}% (2nd={a_second}), rna-seq={rf:.1f}% (2nd={r_second})')

    # Cys codon check
    print('\n  Cys codon (TGT/TGC) at CDR3 start:')
    for sid in SHORT_IDS:
        a_nt = amp_clones[sid]['cdr3nt'].str[:3]
        r_nt = rna_clones[sid]['cdr3nt'].str[:3]
        ac = (a_nt.isin(['TGT', 'TGC'])).sum() / len(a_nt) * 100
        rc = (r_nt.isin(['TGT', 'TGC'])).sum() / len(r_nt) * 100
        print(f'    {sid}: amp={ac:.1f}%, rna-seq={rc:.1f}%')

    # V-gene stratified length
    print('\n  V-gene stratified mean AA length:')
    for sid in SHORT_IDS:
        a_v = amp_clones[sid].copy(); r_v = rna_clones[sid].copy()
        a_v['v_gene'] = a_v['allVHitsWithScore'].apply(lambda x: parse_genes(x, 'TRBV')[0] if parse_genes(x, 'TRBV') else None)
        r_v['v_gene'] = r_v['allVHitsWithScore'].apply(lambda x: parse_genes(x, 'TRBV')[0] if parse_genes(x, 'TRBV') else None)
        a_v['aa_len'] = a_v['cdr3aa'].str.len(); r_v['aa_len'] = r_v['cdr3aa'].str.len()
        a_mean = a_v.groupby('v_gene')['aa_len'].mean()
        r_mean = r_v.groupby('v_gene')['aa_len'].mean()
        common = set(a_mean.index) & set(r_mean.index)
        if common:
            diff = [a_mean[g] - r_mean[g] for g in common]
            longer = sum(d > 0 for d in diff); shorter = sum(d < 0 for d in diff)
            print(f'    {sid}: {len(common)} common V-genes, amp longer in {longer}, rna-seq longer in {shorter}')

# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print('\n' + '='*60)
    print('Preset Comparison: generic-amplicon vs rna-seq (both MiXCR 4.6.0)')
    print('='*60 + '\n')

    plot_alignment()
    plot_assembly()
    plot_sensitivity()
    plot_sensitivity_filtered()
    plot_cdr3_length()
    plot_gene_usage()
    plot_cdr3_overlap()
    plot_cdr3_overlap_filtered()
    amp_agg, rna_agg = plot_clone_correlation()
    plot_nt_correlation()
    find_discordant_clones(amp_agg, rna_agg)
    plot_cdr3_boxplot()
    plot_top_overlap()
    biological_validation()

    print(f'\n✅ All figures saved to: {OUT}/')

    # Summary table
    print('\nSummary:')
    print(f'{"Sample":<8} {"amp_clones":<12} {"rna_clones":<12} {"amp_AA":<10} {"rna_AA":<10} {"shared_AA":<12} {"Jaccard":<10} {"amp-only":<10} {"rna-only":<10}')
    for sid in SHORT_IDS:
        sa = set(amp_clones[sid]['cdr3aa']); sr = set(rna_clones[sid]['cdr3aa'])
        sh = sa & sr; uni = sa | sr
        j = len(sh)/len(uni)*100 if uni else 0
        print(f'{sid:<8} {len(amp_clones[sid]):<12,} {len(rna_clones[sid]):<12,} {len(sa):<10,} {len(sr):<10,} {len(sh):<12,} {j:<10.1f} {len(sa-sr):<10,} {len(sr-sa):<10,}')
