#!/usr/bin/env python3
"""
Diversity and consistency metrics comparing MiXCR generic-amplicon vs rna-seq presets.

Metrics computed per-sample at both CDR3 AA and NT levels:
  Per-preset:   Shannon Entropy, Normalized Shannon (Pielou), Clonality
  Cross-preset: Morisita-Horn Index, Bhattacharyya Coefficient

Quality filters: CDR3 AA length 5-45 aa, readCount > 1
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BASE = "/Users/edwardchan/Desktop/TCR"
OUT  = os.path.join(BASE, "comparison")
os.makedirs(OUT, exist_ok=True)

SAMPLES = [
    "S044_SZ20250508032WHB-2_gdna_genome_2537265",
    "S045_SZ20250522049WHB-1_gdna_genome_2537266",
    "S046_SZ20250522052WHB-9_gdna_genome_2537267",
    "S047_SZ20250522068WHB-0_gdna_genome_2537268",
    "S048_SZ20250621039WHB-2_gdna_genome_2537269",
]
SHORT_IDS = ["S044", "S045", "S046", "S047", "S048"]


def load_clonotypes(path):
    df = pd.read_csv(path, sep='\t')
    df = df.rename(columns={
        'nFeatureSequences(CDR3)': 'cdr3nt',
        'aaFeatureSequences(CDR3)': 'cdr3aa',
    })
    df = df[df['cdr3nt'].notna() & (df['cdr3nt'] != '')]
    df = df[df['cdr3aa'].notna() & (df['cdr3aa'] != '')]
    return df


def filter_clonotypes(df):
    df = df.copy()
    df['aa_len'] = df['cdr3aa'].str.len()
    df = df[(df['aa_len'] >= 5) & (df['aa_len'] <= 45)]
    df = df[df['readCount'] > 1]
    return df


def aggregate_aa(df):
    """Aggregate readCount/readFraction by CDR3 AA, merging synonymous NT variants."""
    agg = df.groupby('cdr3aa').agg(
        readCount=('readCount', 'sum'),
        readFraction=('readFraction', 'sum'),
    ).reset_index()
    return agg


def shannon_entropy(fractions):
    """H = -Σ p_i * ln(p_i)"""
    p = fractions.values
    p = p[p > 0]
    if len(p) == 0:
        return 0.0
    return float(-np.sum(p * np.log(p)))


def normalized_shannon(entropy, n_clones):
    """Pielou's evenness: J' = H / ln(S)"""
    if n_clones <= 1:
        return 1.0
    return entropy / np.log(n_clones)


def clonality(j_prime):
    """1 - J'"""
    return 1.0 - j_prime


def morisita_horn(p_amp, p_rna):
    """C_H = 2 * Σ(p_i * q_i) / (Σp_i² + Σq_i²)

    p_amp and p_rna are pd.Series indexed by clone ID, containing readFraction values.
    Only shared clones contribute to numerator; denominator includes all clones.
    """
    shared = p_amp.index.intersection(p_rna.index)
    num = 2.0 * np.sum(p_amp[shared] * p_rna[shared])
    denom = np.sum(p_amp.values ** 2) + np.sum(p_rna.values ** 2)
    if denom == 0:
        return 0.0
    return float(num / denom)


def bhattacharyya(p_amp, p_rna):
    """BC = Σ √(p_i * q_i) over the union of clones.

    Clones present in only one preset contribute 0 (since one factor is 0).
    """
    shared = p_amp.index.intersection(p_rna.index)
    bc = np.sum(np.sqrt(p_amp[shared] * p_rna[shared]))
    return float(bc)


def compute_indices(amp_df, rna_df, level_label):
    """Compute all diversity and consistency indices for one sample at one level.

    Returns a dict with all metrics.
    """
    # Build frequency series keyed by clone ID
    if level_label == 'AA':
        amp_agg = aggregate_aa(amp_df)
        rna_agg = aggregate_aa(rna_df)
        amp_freq = amp_agg.set_index('cdr3aa')['readFraction']
        rna_freq = rna_agg.set_index('cdr3aa')['readFraction']
        amp_s = len(amp_agg)
        rna_s = len(rna_agg)
    else:
        amp_freq = amp_df.set_index('cdr3nt')['readFraction']
        rna_freq = rna_df.set_index('cdr3nt')['readFraction']
        amp_s = len(amp_df)
        rna_s = len(rna_df)

    # Per-preset diversity
    amp_h = shannon_entropy(amp_freq)
    rna_h = shannon_entropy(rna_freq)
    amp_j = normalized_shannon(amp_h, amp_s)
    rna_j = normalized_shannon(rna_h, rna_s)

    # Cross-preset consistency
    mh = morisita_horn(amp_freq, rna_freq)
    bc = bhattacharyya(amp_freq, rna_freq)

    return {
        'level': level_label,
        'amp_clones': amp_s,
        'rna_clones': rna_s,
        'amp_shannon_H': round(amp_h, 4),
        'rna_shannon_H': round(rna_h, 4),
        'amp_pielou_J': round(amp_j, 4),
        'rna_pielou_J': round(rna_j, 4),
        'amp_clonality': round(clonality(amp_j), 4),
        'rna_clonality': round(clonality(rna_j), 4),
        'morisita_horn': round(mh, 6),
        'bhattacharyya_BC': round(bc, 4),
    }


def plot_heatmap(results_df):
    """Heatmap of all 5 samples × selected metrics at both levels."""
    metrics = [
        'amp_shannon_H', 'rna_shannon_H',
        'amp_pielou_J', 'rna_pielou_J',
        'morisita_horn', 'bhattacharyya_BC',
    ]
    labels = [
        'amp\nShannon H', 'rna-seq\nShannon H',
        'amp\nPielou J\'', 'rna-seq\nPielou J\'',
        'Morisita-\nHorn', 'Bhattacharyya\nBC',
    ]

    # Pivot: rows = sample+level, columns = metric
    pivot_data = {}
    for _, row in results_df.iterrows():
        key = f"{row['sample']}_{row['level']}"
        pivot_data[key] = [row[m] for m in metrics]

    plot_df = pd.DataFrame(pivot_data, index=labels).T

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.heatmap(plot_df, annot=True, fmt='.3f', cmap='YlOrRd', ax=ax,
                linewidths=0.5, cbar_kws={'label': 'Value'})
    ax.set_title('P11: Diversity & Consistency Metrics\n(generic-amplicon vs rna-seq, filtered)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    fig_dir = os.path.join(OUT, 'figures_preset')
    os.makedirs(fig_dir, exist_ok=True)
    fig.savefig(f'{fig_dir}/P11_diversity_heatmap.png', dpi=150)
    plt.close()
    print(f'  ✓ P11_diversity_heatmap.png')


def main():
    print('\n' + '=' * 60)
    print('Diversity & Consistency Metrics')
    print('Filters: CDR3 AA length 5-45, readCount > 1')
    print('=' * 60 + '\n')

    results = []

    for sample, sid in zip(SAMPLES, SHORT_IDS):
        amp_raw = load_clonotypes(
            f'{BASE}/test_v3/{sample}/Map_Clone_Analysis/{sample}.clonotypes.TRB.raw.txt')
        rna_raw = load_clonotypes(
            f'{BASE}/test_v3_rnaseq/{sample}/Map_Clone_Analysis/{sample}.clonotypes.TRB.raw.txt')

        amp_filt = filter_clonotypes(amp_raw)
        rna_filt = filter_clonotypes(rna_raw)

        # AA level
        aa_results = compute_indices(amp_filt, rna_filt, 'AA')
        aa_results['sample'] = sid
        results.append(aa_results)

        # NT level
        nt_results = compute_indices(amp_filt, rna_filt, 'NT')
        nt_results['sample'] = sid
        results.append(nt_results)

        print(f'  {sid} (filtered): amp={len(amp_filt):,} clones, rna-seq={len(rna_filt):,} clones')
        print(f'    AA: H_amp={aa_results["amp_shannon_H"]:.4f}, H_rna={aa_results["rna_shannon_H"]:.4f}, '
              f'MH={aa_results["morisita_horn"]:.6f}, BC={aa_results["bhattacharyya_BC"]:.4f}')
        print(f'    NT: H_amp={nt_results["amp_shannon_H"]:.4f}, H_rna={nt_results["rna_shannon_H"]:.4f}, '
              f'MH={nt_results["morisita_horn"]:.6f}, BC={nt_results["bhattacharyya_BC"]:.4f}')

    results_df = pd.DataFrame(results)
    cols = ['sample', 'level', 'amp_clones', 'rna_clones',
            'amp_shannon_H', 'rna_shannon_H',
            'amp_pielou_J', 'rna_pielou_J',
            'amp_clonality', 'rna_clonality',
            'morisita_horn', 'bhattacharyya_BC']
    results_df = results_df[cols]

    tsv_path = os.path.join(OUT, 'diversity_metrics.tsv')
    results_df.to_csv(tsv_path, sep='\t', index=False)
    print(f'\n  ✓ Results saved to {tsv_path}')

    plot_heatmap(results_df)
    print('\n✅ Done.')


if __name__ == '__main__':
    main()
