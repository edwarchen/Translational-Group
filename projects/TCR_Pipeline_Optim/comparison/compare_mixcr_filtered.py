#!/usr/bin/env python3
"""
MiXCR preset comparison — FILTERED ANALYSIS
Filters: readCount > 5, CDR3 AA length > 5 and < 45
Compares v1 (3.0.4 / rna-seq) vs v3 (4.6.0 / generic-amplicon --dna) after filtering.
"""

import os, re
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# ─── Configuration ───────────────────────────────────────────────────────────

BASE = "/Users/edwardchan/Desktop/TCR"
OUT  = os.path.join(BASE, "comparison", "figures_filtered")
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

# ─── Helper functions ───────────────────────────────────────────────────────

def load_clonotypes(path, version='v1'):
    df = pd.read_csv(path, sep='\t')
    if version == 'v1':
        df = df.rename(columns={
            'cloneCount': 'readCount', 'cloneFraction': 'readFraction',
            'nSeqCDR3': 'cdr3nt', 'aaSeqCDR3': 'cdr3aa',
        })
    else:
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
        if m:
            genes.add(m.group(1))
    return sorted(genes)

# ─── Load raw data ──────────────────────────────────────────────────────────

print('Loading raw data...')
v1_raw, v3_raw = {}, {}
for s in SAMPLES:
    sid = SHORT_MAP[s]
    v1_raw[sid] = load_clonotypes(
        f'{BASE}/test_v1/{s}/Map_Clone_Analysis/{s}.clonotypes.TRB.raw.txt', 'v1')
    v3_raw[sid] = load_clonotypes(
        f'{BASE}/test_v3/{s}/Map_Clone_Analysis/{s}.clonotypes.TRB.raw.txt', 'v3')
    print(f'  {sid}: v1 raw={len(v1_raw[sid]):,} clones, v3 raw={len(v3_raw[sid]):,} clones')

# ─── Apply filters ──────────────────────────────────────────────────────────

print('\nApplying filters: readCount > 5, CDR3 AA length > 5 and < 45...')
v1_f, v3_f = {}, {}

for sid in SHORT_IDS:
    v1d = v1_raw[sid].copy()
    v3d = v3_raw[sid].copy()

    v1d['aa_len'] = v1d['cdr3aa'].str.len()
    v3d['aa_len'] = v3d['cdr3aa'].str.len()

    v1_pre = len(v1d)
    v3_pre = len(v3d)

    v1d = v1d[(v1d['readCount'] > 5) & (v1d['aa_len'] > 5) & (v1d['aa_len'] < 45)]
    v3d = v3d[(v3d['readCount'] > 5) & (v3d['aa_len'] > 5) & (v3d['aa_len'] < 45)]

    v1_f[sid] = v1d
    v3_f[sid] = v3d

    print(f'  {sid}: v1 {v1_pre:,} → {len(v1d):,} ({len(v1d)/v1_pre*100:.1f}%), '
          f'v3 {v3_pre:,} → {len(v3d):,} ({len(v3d)/v3_pre*100:.1f}%)')

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE F1: Before/After Filtering — Clone Count Change
# ═════════════════════════════════════════════════════════════════════════════

def plot_filter_effect():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    x = np.arange(len(SHORT_IDS))
    w = 0.2
    C1, C2 = '#4C72B0', '#DD8452'

    # A: Clone count before/after
    ax = axes[0]
    for i, sid in enumerate(SHORT_IDS):
        raw1, raw3 = len(v1_raw[sid]), len(v3_raw[sid])
        fil1, fil3 = len(v1_f[sid]), len(v3_f[sid])
        bars = ax.bar(i - w, raw1, w*0.8, color=C1, alpha=0.3, label='v1 raw' if i == 0 else '')
        ax.bar(i - w, fil1, w*0.8, color=C1, alpha=1.0, label='v1 filtered' if i == 0 else '')
        ax.bar(i + w, raw3, w*0.8, color=C2, alpha=0.3, label='v3 raw' if i == 0 else '')
        ax.bar(i + w, fil3, w*0.8, color=C2, alpha=1.0, label='v3 filtered' if i == 0 else '')
    ax.set_ylabel('Clone Count')
    ax.set_title('A. Clone Count: Raw (light) vs Filtered (solid)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)
    ax.legend(fontsize=8)

    # B: Survival rate %
    ax = axes[1]
    s1 = [len(v1_f[s]) / len(v1_raw[s]) * 100 for s in SHORT_IDS]
    s3 = [len(v3_f[s]) / len(v3_raw[s]) * 100 for s in SHORT_IDS]
    ax.bar(x - w/2, s1, w, color=C1, label='v1')
    ax.bar(x + w/2, s3, w, color=C2, label='v3')
    ax.set_ylabel('Clones Surviving Filter (%)')
    ax.set_title('B. Filter Survival Rate (%)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)
    ax.legend(fontsize=9)
    for i in range(len(SHORT_IDS)):
        ax.text(i - w/2, s1[i] + 0.5, f'{s1[i]:.1f}%', ha='center', fontsize=8)
        ax.text(i + w/2, s3[i] + 0.5, f'{s3[i]:.1f}%', ha='center', fontsize=8)

    fig.suptitle('F1: Effect of Filtering (readCount>5, AA 5-45)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/F1_filter_effect.png')
    plt.close()
    print('  ✓ F1_filter_effect.png')

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE F2: CDR3 AA Overlap — Filtered
# ═════════════════════════════════════════════════════════════════════════════

def plot_cdr3_overlap_filtered():
    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(SHORT_IDS))
    w = 0.25

    v1_only, v3_only, shared = [], [], []
    v1_pct, v3_pct, sh_pct = [], [], []

    for sid in SHORT_IDS:
        s1 = set(v1_f[sid]['cdr3aa'])
        s3 = set(v3_f[sid]['cdr3aa'])
        sh = s1 & s3
        total = len(s1 | s3)
        v1_only.append(len(s1 - s3))
        v3_only.append(len(s3 - s1))
        shared.append(len(sh))
        v1_pct.append(len(s1 - s3) / total * 100 if total else 0)
        v3_pct.append(len(s3 - s1) / total * 100 if total else 0)
        sh_pct.append(len(sh) / total * 100 if total else 0)

    ax.bar(x, v1_only, w, label='v1-only CDR3 AA', color='#4C72B0')
    ax.bar(x, shared, w, bottom=v1_only, label='Shared CDR3 AA', color='#55A868')
    bot = [a + b for a, b in zip(v1_only, shared)]
    ax.bar(x, v3_only, w, bottom=bot, label='v3-only CDR3 AA', color='#DD8452')

    for i in range(len(SHORT_IDS)):
        ax.text(i, v1_only[i]/2, f'{v1_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(i, v1_only[i] + shared[i]/2, f'{sh_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')
        ax.text(i, v1_only[i] + shared[i] + v3_only[i]/2, f'{v3_pct[i]:.1f}%', ha='center', va='center', fontsize=9, fontweight='bold')

    ax.set_ylabel('CDR3 AA Sequence Count')
    ax.set_title('CDR3 AA Overlap After Filtering (readCount>5, AA 5-45)')
    ax.set_xticks(x); ax.set_xticklabels(SHORT_IDS)
    ax.legend(fontsize=10)

    fig.suptitle('F2: CDR3 Sequence Overlap — Filtered', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/F2_cdr3_overlap_filtered.png')
    plt.close()
    print('  ✓ F2_cdr3_overlap_filtered.png')

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE F3: Shared Clone Frequency Correlation — Filtered + Stratified
# ═════════════════════════════════════════════════════════════════════════════

def plot_clone_correlation_filtered():
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for i, sid in enumerate(SHORT_IDS):
        ax = axes[i]
        v1d = v1_f[sid][['cdr3aa', 'readFraction', 'readCount']].copy()
        v3d = v3_f[sid][['cdr3aa', 'readFraction', 'readCount']].copy()
        merged = v1d.merge(v3d, on='cdr3aa', suffixes=('_v1', '_v3'))

        n_shared = len(merged)
        # Count v1-only and v3-only in filtered sets
        s1 = set(v1d['cdr3aa']); s3 = set(v3d['cdr3aa'])
        n_v1_only = len(s1 - s3)
        n_v3_only = len(s3 - s1)

        if len(merged) > 0:
            ax.scatter(merged['readFraction_v1'] * 100, merged['readFraction_v3'] * 100,
                      alpha=0.3, s=8, edgecolors='none')
            ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlabel('v1 Frequency (%)'); ax.set_ylabel('v3 Frequency (%)')
        rho = merged['readFraction_v1'].corr(merged['readFraction_v3'], method='spearman') if len(merged) > 1 else 0
        ax.set_title(f'{sid}\nShared: {n_shared:,} | v1-only: {n_v1_only:,} | v3-only: {n_v3_only:,}\nSpearman ρ = {rho:.4f}')
        ax.plot([1e-6, 100], [1e-6, 100], 'k--', alpha=0.3, linewidth=1)

    axes[5].set_visible(False)
    fig.suptitle('F3: Clone Frequency Correlation — Filtered (readCount>5, AA 5-45)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/F3_clone_correlation_filtered.png')
    plt.close()
    print('  ✓ F3_clone_correlation_filtered.png')

    # ── Stratified correlation table ──
    print('\n  Stratified frequency correlation (filtered shared clones):')
    for sid in SHORT_IDS:
        m = v1_f[sid][['cdr3aa','readFraction','readCount']].merge(
            v3_f[sid][['cdr3aa','readFraction','readCount']], on='cdr3aa', suffixes=('_v1','_v3'))
        m['mean_reads'] = (m['readCount_v1'] + m['readCount_v3']) / 2
        bins = [(6, 10), (11, 50), (51, 200), (201, 1000), (1001, 999999)]
        print(f'    {sid} (shared={len(m):,}):')
        for lo, hi in bins:
            sub = m[(m['mean_reads'] >= lo) & (m['mean_reads'] <= hi)]
            if len(sub) < 3: continue
            rho = sub['readFraction_v1'].corr(sub['readFraction_v3'], method='spearman')
            r   = sub['readFraction_v1'].corr(sub['readFraction_v3'], method='pearson')
            ratio = sub['readFraction_v1'] / sub['readFraction_v3']
            gt2x = ((ratio > 2) | (ratio < 0.5)).sum()
            print(f'      {lo}-{hi:>5} reads:  N={len(sub):,}  ρ={rho:+.3f}  r={r:+.3f}  >2x={gt2x/len(sub)*100:.1f}%')

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE F4: CDR3 Length Distribution — Filtered
# ═════════════════════════════════════════════════════════════════════════════

def plot_cdr3_length_filtered():
    fig, axes = plt.subplots(2, 1, figsize=(16, 10))

    # A: NT length
    ax = axes[0]
    for i, sid in enumerate(SHORT_IDS):
        l1 = v1_f[sid]['cdr3nt'].str.len().dropna()
        l3 = v3_f[sid]['cdr3nt'].str.len().dropna()
        bins = np.arange(0, 150, 3)
        ax.hist(l1, bins=bins, alpha=0.45, color=SAMPLE_COLORS[i], label=f'{sid} v1',
                density=True, histtype='stepfilled', linewidth=1.5, edgecolor=SAMPLE_COLORS[i])
        ax.hist(l3, bins=bins, alpha=0.25, color=SAMPLE_COLORS[i], label=f'{sid} v3',
                density=True, histtype='step', linewidth=2.5, linestyle='--',
                edgecolor=SAMPLE_COLORS[i])
    ax.set_xlabel('CDR3 Nucleotide Length (bp)')
    ax.set_ylabel('Density')
    ax.set_title('A. CDR3 NT Length Distribution — Filtered (solid=v1, dashed=v3)')
    ax.set_xlim(15, 150)
    ax.legend(fontsize=8, ncol=2)

    # B: AA length
    ax = axes[1]
    for i, sid in enumerate(SHORT_IDS):
        l1 = v1_f[sid]['cdr3aa'].str.len().dropna()
        l3 = v3_f[sid]['cdr3aa'].str.len().dropna()
        bins = np.arange(5, 45, 1)
        ax.hist(l1, bins=bins, alpha=0.45, color=SAMPLE_COLORS[i], label=f'{sid} v1',
                density=True, histtype='stepfilled', linewidth=1.5, edgecolor=SAMPLE_COLORS[i])
        ax.hist(l3, bins=bins, alpha=0.25, color=SAMPLE_COLORS[i], label=f'{sid} v3',
                density=True, histtype='step', linewidth=2.5, linestyle='--',
                edgecolor=SAMPLE_COLORS[i])
    ax.set_xlabel('CDR3 Amino Acid Length')
    ax.set_ylabel('Density')
    ax.set_title('B. CDR3 AA Length Distribution — Filtered (solid=v1, dashed=v3)')
    ax.set_xlim(5, 45)
    ax.legend(fontsize=8, ncol=2)

    fig.suptitle('F4: CDR3 Length Distribution — Filtered (readCount>5, AA 5-45)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/F4_cdr3_length_filtered.png')
    plt.close()
    print('  ✓ F4_cdr3_length_filtered.png')

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE F5: V/J Gene Usage — Filtered
# ═════════════════════════════════════════════════════════════════════════════

def plot_gene_usage_filtered():
    fig, axes = plt.subplots(2, 1, figsize=(16, 12))

    for gene_type, ax in [('V', axes[0]), ('J', axes[1])]:
        prefix = 'TRB' + gene_type
        hit_col = f'all{gene_type}HitsWithScore'

        v1_counts = Counter()
        v3_counts = Counter()
        for sid in SHORT_IDS:
            for _, row in v1_f[sid].iterrows():
                for g in parse_genes(row.get(hit_col, ''), prefix):
                    v1_counts[g] += row['readCount']
            for _, row in v3_f[sid].iterrows():
                for g in parse_genes(row.get(hit_col, ''), prefix):
                    v3_counts[g] += row['readCount']

        v1_total = sum(v1_counts.values())
        v3_total = sum(v3_counts.values())
        all_genes = sorted(set(v1_counts) | set(v3_counts))

        f1 = [v1_counts.get(g, 0) / v1_total * 100 for g in all_genes]
        f3 = [v3_counts.get(g, 0) / v3_total * 100 for g in all_genes]

        ax.scatter(f1, f3, alpha=0.6, s=50, edgecolors='grey', linewidth=0.5)
        mx = max(max(f1), max(f3)) * 1.05
        ax.plot([0, mx], [0, mx], 'k--', alpha=0.3, linewidth=1)

        diffs = [abs(f1[i] - f3[i]) for i in range(len(all_genes))]
        for idx in np.argsort(diffs)[-10:]:
            ax.annotate(all_genes[idx], (f1[idx], f3[idx]),
                       fontsize=7, alpha=0.8, xytext=(5, 5), textcoords='offset points')

        rho = pd.Series(f1).corr(pd.Series(f3), method='spearman')
        ax.set_xlabel(f'v1 {gene_type} Gene Frequency (%)')
        ax.set_ylabel(f'v3 {gene_type} Gene Frequency (%)')
        ax.set_title(f'{gene_type} Gene Usage — Filtered (Spearman ρ = {rho:.4f})')

    fig.suptitle('F5: V/J Gene Usage — Filtered (readCount>5, AA 5-45)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/F5_gene_usage_filtered.png')
    plt.close()
    print('  ✓ F5_gene_usage_filtered.png')

# ═════════════════════════════════════════════════════════════════════════════
# FIGURE F6: Top-N Clone Overlap — Filtered
# ═════════════════════════════════════════════════════════════════════════════

def plot_top_overlap_filtered():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    top_ns = [10, 50, 100, 500, 1000]
    x = np.arange(len(top_ns))

    ax = axes[0]
    for i, sid in enumerate(SHORT_IDS):
        overlaps = []
        for n in top_ns:
            t1 = set(v1_f[sid].nlargest(min(n, len(v1_f[sid])), 'readCount')['cdr3aa'])
            t3 = set(v3_f[sid].nlargest(min(n, len(v3_f[sid])), 'readCount')['cdr3aa'])
            overlaps.append(len(t1 & t3) / n * 100 if n > 0 else 0)
        ax.plot(x, overlaps, 'o-', color=SAMPLE_COLORS[i], label=sid, linewidth=2, markersize=8)
    ax.set_xticks(x); ax.set_xticklabels([f'Top {n}' for n in top_ns])
    ax.set_ylabel('Overlap (%)')
    ax.set_title('A. Top-N Clone Overlap — Filtered')
    ax.legend(); ax.set_ylim(0, 105)

    ax = axes[1]
    for i, sid in enumerate(SHORT_IDS):
        fracs = []
        for n in top_ns:
            t1 = set(v1_f[sid].nlargest(min(n, len(v1_f[sid])), 'readCount')['cdr3aa'])
            t3 = set(v3_f[sid].nlargest(min(n, len(v3_f[sid])), 'readCount')['cdr3aa'])
            sh = t1 & t3
            tot = v1_f[sid]['readCount'].sum()
            sr = v1_f[sid][v1_f[sid]['cdr3aa'].isin(sh)]['readCount'].sum()
            fracs.append(sr / tot * 100 if tot > 0 else 0)
        ax.plot(x, fracs, 'o-', color=SAMPLE_COLORS[i], label=sid, linewidth=2, markersize=8)
    ax.set_xticks(x); ax.set_xticklabels([f'Top {n}' for n in top_ns])
    ax.set_ylabel('v1 Read Fraction in Shared Clones (%)')
    ax.set_title('B. Read Fraction Captured by Shared Top-N — Filtered')
    ax.legend()

    fig.suptitle('F6: Top Clones Concordance — Filtered (readCount>5, AA 5-45)', fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig(f'{OUT}/F6_top_overlap_filtered.png')
    plt.close()
    print('  ✓ F6_top_overlap_filtered.png')

# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print('\n' + '='*60)
    print('Filtered Analysis: readCount > 5, CDR3 AA length 5-45')
    print('='*60 + '\n')

    plot_filter_effect()
    plot_cdr3_overlap_filtered()
    plot_clone_correlation_filtered()
    plot_cdr3_length_filtered()
    plot_gene_usage_filtered()
    plot_top_overlap_filtered()

    print(f'\n✅ All filtered figures saved to: {OUT}/')

    # ── Print summary table ──
    print('\nSummary Table (filtered):')
    print(f'{"Sample":<8} {"v1_clones":<12} {"v3_clones":<12} {"v1_AA":<10} {"v3_AA":<10} {"shared_AA":<12} {"Jaccard":<10} {"v1-only":<10} {"v3-only":<10}')
    for sid in SHORT_IDS:
        s1 = set(v1_f[sid]['cdr3aa']); s3 = set(v3_f[sid]['cdr3aa'])
        sh = s1 & s3; uni = s1 | s3
        j = len(sh)/len(uni)*100 if uni else 0
        print(f'{sid:<8} {len(v1_f[sid]):<12,} {len(v3_f[sid]):<12,} {len(s1):<10,} {len(s3):<10,} {len(sh):<12,} {j:<10.1f} {len(s1-s3):<10,} {len(s3-s1):<10,}')
