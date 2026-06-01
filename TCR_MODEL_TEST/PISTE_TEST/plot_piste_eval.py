#!/usr/bin/env python3
"""
Nature-style ROC and PR curves for PISTE evaluation on external VDJdb test sets.
"""

import os, sys, random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.utils.data as Data
import torch.nn.functional as F
from sklearn.metrics import roc_curve, precision_recall_curve, auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# ── Paths ──────────────────────────────────────────────
PISTE_CODE = '/home/chenya/TCR_MODEL_TEST/PISTE/code'
VDJDB_DIR = '/home/chenya/TCR_MODEL_TEST/data'
HLA_SEQ_FILE = '/home/chenya/TCR_MODEL_TEST/PISTE/data/raw_data/common_hla_sequence.csv'
OUTPUT_DIR = '/home/chenya/TCR_MODEL_TEST/data/piste_eval'
os.makedirs(OUTPUT_DIR, exist_ok=True)

sys.path.insert(0, PISTE_CODE)
from Model.PISTE import Transformer

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)

HLA_EXPAND = {
    'HLA-A01': 'HLA-A01:01', 'HLA-A02': 'HLA-A02:01', 'HLA-A03': 'HLA-A03:01',
    'HLA-A11': 'HLA-A11:01', 'HLA-B07': 'HLA-B07:02', 'HLA-B08': 'HLA-B08:01',
    'HLA-B27': 'HLA-B27:05', 'HLA-B35': 'HLA-B35:01', 'HLA-B42': 'HLA-B42:01',
    'HLA-B57': 'HLA-B57:01', 'HLA-B58': 'HLA-B58:01',
}

# ── Nature style setup ─────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 8,
    'axes.titlesize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 6.5,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'lines.linewidth': 1.0,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
})

# Nature palette — colorblind-friendly, muted
COLORS = {
    'random': '#377EB8',   # blue
    'unipep': '#4DAF4A',   # green
    'reftcr': '#984EA3',   # purple
}
MODEL_LABELS = {
    'random': 'PISTE-random',
    'unipep': 'PISTE-unipep',
    'reftcr': 'PISTE-reftcr',
}


def load_vdjdb(filepath):
    df = pd.read_csv(filepath, sep='\t')
    df = df.rename(columns={'cdr3.beta': 'CDR3', 'peptide': 'MT_pep'})
    df['HLA_type'] = df['hla'].str.replace('*', '', regex=False)
    df['HLA_type'] = df['HLA_type'].replace(HLA_EXPAND)
    return df[['CDR3', 'MT_pep', 'HLA_type']]


def generate_negatives(pos_df, neg_ratio=10):
    cdrs = pos_df['CDR3'].values
    peps = pos_df['MT_pep'].values
    hlas = pos_df['HLA_type'].values
    pos_set = set(zip(cdrs, peps, hlas))
    negs = []
    for _ in range(len(pos_df) * neg_ratio):
        while True:
            t = (random.choice(cdrs), random.choice(peps), random.choice(hlas))
            if t not in pos_set:
                negs.append(t)
                pos_set.add(t)
                break
    return pd.DataFrame(negs, columns=['CDR3', 'MT_pep', 'HLA_type'])


def run_prediction(pos_df, model_name, neg_ratio=10):
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    neg_df = generate_negatives(pos_df, neg_ratio)
    pos_df, neg_df = pos_df.copy(), neg_df.copy()
    pos_df['Label'], neg_df['Label'] = 1, 0
    df = pd.concat([pos_df, neg_df], ignore_index=True).sample(frac=1, random_state=42)

    hla_seq = pd.read_csv(HLA_SEQ_FILE)
    df = pd.merge(df, hla_seq, on='HLA_type', how='inner')

    vocab = {'C':1,'W':2,'V':3,'A':4,'H':5,'T':6,'E':7,'K':8,'N':9,'P':10,
             'I':11,'L':12,'S':13,'D':14,'G':15,'Q':16,'R':17,'Y':18,'F':19,'M':20,'-':0}
    pep_max_len, hla_max_len, tcr_max_len = 11, 34, 30
    pep_arr, hla_arr, tcr_arr, lbl_arr = [], [], [], []
    for _, r in df.iterrows():
        pep_arr.append([vocab.get(aa,0) for aa in str(r['MT_pep']).ljust(pep_max_len,'-')])
        hla_arr.append([vocab.get(aa,0) for aa in str(r['HLA_sequence']).ljust(hla_max_len,'-')])
        tcr_arr.append([vocab.get(aa,0) for aa in str(r['CDR3']).ljust(tcr_max_len,'-')])
        lbl_arr.append(int(r['Label']))

    class DS(Data.Dataset):
        def __init__(self, p, h, t, l):
            super().__init__()
            self.p, self.h, self.t, self.l = p, h, t, l
        def __len__(self):
            return self.p.shape[0]
        def __getitem__(self, i):
            return self.p[i], self.h[i], self.t[i], self.l[i]

    loader = Data.DataLoader(
        DS(torch.LongTensor(np.array(pep_arr, dtype=np.int64)),
           torch.LongTensor(np.array(hla_arr, dtype=np.int64)),
           torch.LongTensor(np.array(tcr_arr, dtype=np.int64)),
           torch.LongTensor(lbl_arr)),
        batch_size=1024, shuffle=False, num_workers=0)

    window_map = {'random': '3', 'unipep': '1', 'reftcr': '2'}
    model = Transformer(
        device=device, vocab_size=21, d_model=64, e_layers=3, d=64, n_heads=9,
        sigma=1, window_threshold=window_map.get(model_name,'3'), d_ff=512,
        interact_layers=1, tgt_len=75, hla_max_len=34, d_layers=1).to(device)
    ckpt = os.path.join(PISTE_CODE, 'checkpoints', model_name, 'exp0.pkl')
    model.load_state_dict(torch.load(ckpt, map_location=device))
    model.eval()

    all_lbls, all_probs = [], []
    with torch.no_grad():
        for bp, bh, bt, bl in loader:
            bp, bh, bt = bp.to(device), bh.to(device), bt.to(device)
            out, _, _ = model(bp, bh, bt)
            all_lbls.extend(bl.numpy())
            all_probs.extend(F.softmax(out, dim=1)[:, 1].cpu().numpy())

    return np.array(all_lbls), np.array(all_probs)


def plot_curves(results, save_path):
    """Create a 2×2 panel figure: ROC (top) and PR (bottom) for both datasets."""
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.8))

    dataset_labels = {
        'score_gt_0': 'Independent VDJdb (score > 0)',
        'score_ge_2': 'Independent VDJdb (score ≥ 2)',
    }

    # ── ROC curves (top row) ──
    for col, dname in enumerate(['score_gt_0', 'score_ge_2']):
        ax = axes[0, col]
        for mname in ['random', 'unipep', 'reftcr']:
            y_true, y_prob = results[dname][mname]
            fpr, tpr, _ = roc_curve(y_true, y_prob)
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, color=COLORS[mname], label=f'{MODEL_LABELS[mname]} (AUC={roc_auc:.3f})',
                    rasterized=False)

        ax.plot([0, 1], [0, 1], '--', color='grey', linewidth=0.5, alpha=0.7)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.legend(loc='lower right', frameon=False, handlelength=1.2,
                  handletextpad=0.5, borderpad=0.3)
        ax.set_title(dataset_labels[dname], fontweight='normal', pad=4)
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_major_locator(MultipleLocator(0.2))

    # ── PR curves (bottom row) ──
    for col, dname in enumerate(['score_gt_0', 'score_ge_2']):
        ax = axes[1, col]
        for mname in ['random', 'unipep', 'reftcr']:
            y_true, y_prob = results[dname][mname]
            prec, rec, _ = precision_recall_curve(y_true, y_prob)
            pr_auc = auc(rec, prec)
            ax.plot(rec, prec, color=COLORS[mname], label=f'{MODEL_LABELS[mname]} (AUC={pr_auc:.3f})')

        # Baseline: random classifier
        n_pos = sum(results[dname]['random'][0])
        n_total = len(results[dname]['random'][0])
        baseline = n_pos / n_total
        ax.axhline(y=baseline, color='grey', linestyle='--', linewidth=0.5, alpha=0.7)

        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')
        ax.legend(loc='upper right', frameon=False, handlelength=1.2,
                  handletextpad=0.5, borderpad=0.3)
        ax.set_title(dataset_labels[dname], fontweight='normal', pad=4)
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.yaxis.set_major_locator(MultipleLocator(0.2))

    # ── Panel labels ──
    for idx, (ax, label) in enumerate(zip(axes.flat, ['a', 'b', 'c', 'd'])):
        ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=9,
                fontweight='bold', va='top', ha='left')

    plt.tight_layout(pad=1.2, h_pad=1.5, w_pad=1.5)
    fig.savefig(save_path, dpi=300, bbox_inches='tight', pad_inches=0.08,
                facecolor='white', edgecolor='none')
    fig.savefig(save_path.replace('.pdf', '.png'), dpi=300, bbox_inches='tight',
                pad_inches=0.08, facecolor='white', edgecolor='none')
    print(f'Saved: {save_path}')
    plt.close()


def main():
    datasets = {
        'score_gt_0': os.path.join(VDJDB_DIR, 'vdjdb_filtered_score_gt_0_clean.tsv'),
        'score_ge_2': os.path.join(VDJDB_DIR, 'vdjdb_filtered_score_ge_2_clean.tsv'),
    }

    results = {}
    for dname, fpath in datasets.items():
        print(f'\n{"="*60}')
        print(f'Dataset: {dname}')
        pos_df = load_vdjdb(fpath)
        hla_db = set(pd.read_csv(HLA_SEQ_FILE)['HLA_type'])
        unmatched = set(pos_df['HLA_type']) - hla_db
        if unmatched:
            pos_df = pos_df[~pos_df['HLA_type'].isin(unmatched)]
        print(f'  {len(pos_df)} positive triples')

        results[dname] = {}
        for mname in ['random', 'unipep', 'reftcr']:
            print(f'  PISTE-{mname}...')
            y_true, y_prob = run_prediction(pos_df, mname, neg_ratio=10)
            results[dname][mname] = (y_true, y_prob)
            roc_auc = auc(*roc_curve(y_true, y_prob)[:2])
            prec, rec, _ = precision_recall_curve(y_true, y_prob)
            aupr = auc(rec, prec)
            print(f'    AUC={roc_auc:.4f}  AUPR={aupr:.4f}')

            # Save raw predictions
            pred_file = os.path.join(OUTPUT_DIR, f'{dname}_{mname}_preds.npz')
            np.savez(pred_file, y_true=y_true, y_prob=y_prob)

    # Plot
    print(f"\n{'='*60}")
    print('Generating Nature-style figures...')
    plot_curves(results, os.path.join(OUTPUT_DIR, 'piste_evaluation_curves.pdf'))


if __name__ == '__main__':
    main()
