#!/usr/bin/env python3
"""
Evaluate PISTE on external VDJdb test sets.
Run with: conda activate piste && python3 run_piste_eval.py
"""

import os, sys, random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.utils.data as Data
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

# Paths
PISTE_CODE = '/home/chenya/TCR_MODEL_TEST/PISTE/code'
VDJDB_DIR = '/home/chenya/TCR_MODEL_TEST/data'
HLA_SEQ_FILE = '/home/chenya/TCR_MODEL_TEST/PISTE/data/raw_data/common_hla_sequence.csv'
OUTPUT_DIR = '/home/chenya/TCR_MODEL_TEST/data/piste_eval'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Import PISTE Transformer
sys.path.insert(0, PISTE_CODE)
from Model.PISTE import Transformer

random.seed(42)
np.random.seed(42)
torch.manual_seed(42)

# HLA mapping: 2-field -> most common 4-field subtype in PISTE DB
HLA_EXPAND = {
    'HLA-A01': 'HLA-A01:01',
    'HLA-A02': 'HLA-A02:01',
    'HLA-A03': 'HLA-A03:01',
    'HLA-A11': 'HLA-A11:01',
    'HLA-B07': 'HLA-B07:02',
    'HLA-B08': 'HLA-B08:01',
    'HLA-B27': 'HLA-B27:05',
    'HLA-B35': 'HLA-B35:01',
    'HLA-B42': 'HLA-B42:01',
    'HLA-B57': 'HLA-B57:01',
    'HLA-B58': 'HLA-B58:01',
}


def load_vdjdb(filepath):
    """Load filtered VDJdb, convert to PISTE format."""
    df = pd.read_csv(filepath, sep='\t')
    df = df.rename(columns={'cdr3.beta': 'CDR3', 'peptide': 'MT_pep'})
    # Convert HLA: HLA-A*02:01 -> HLA-A02:01
    df['HLA_type'] = df['hla'].str.replace('*', '', regex=False)
    # Expand 2-field HLAs to 4-field
    df['HLA_type'] = df['HLA_type'].replace(HLA_EXPAND)
    return df[['CDR3', 'MT_pep', 'HLA_type']]


def generate_negatives(pos_df, neg_ratio=10):
    """Random-shuffle negative sampling."""
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
    """Run PISTE on a single dataset."""
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    # Prepare data with negatives
    neg_df = generate_negatives(pos_df, neg_ratio)
    pos_df_labeled = pos_df.copy()
    pos_df_labeled['Label'] = 1
    neg_df_labeled = neg_df.copy()
    neg_df_labeled['Label'] = 0
    df = pd.concat([pos_df_labeled, neg_df_labeled], ignore_index=True)
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)

    # Merge HLA sequences
    hla_seq = pd.read_csv(HLA_SEQ_FILE)
    n_before = len(df)
    df = pd.merge(df, hla_seq, on='HLA_type', how='inner')
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(f"    Dropped {n_dropped} triples with unmapped HLA (before merge)")

    # Make tokenized inputs
    vocab = {'C':1,'W':2,'V':3,'A':4,'H':5,'T':6,'E':7,'K':8,'N':9,'P':10,
             'I':11,'L':12,'S':13,'D':14,'G':15,'Q':16,'R':17,'Y':18,'F':19,'M':20,'-':0}
    pep_max_len, hla_max_len, tcr_max_len = 11, 34, 30
    tgt_len = pep_max_len + hla_max_len + tcr_max_len

    pep_arr, hla_arr, tcr_arr, lbl_arr = [], [], [], []
    for _, r in df.iterrows():
        pep_arr.append([vocab.get(aa,0) for aa in str(r['MT_pep']).ljust(pep_max_len,'-')])
        hla_arr.append([vocab.get(aa,0) for aa in str(r['HLA_sequence']).ljust(hla_max_len,'-')])
        tcr_arr.append([vocab.get(aa,0) for aa in str(r['CDR3']).ljust(tcr_max_len,'-')])
        lbl_arr.append(int(r['Label']))

    pep_t = torch.LongTensor(np.array(pep_arr, dtype=np.int64))
    hla_t = torch.LongTensor(np.array(hla_arr, dtype=np.int64))
    tcr_t = torch.LongTensor(np.array(tcr_arr, dtype=np.int64))
    lbl_t = torch.LongTensor(lbl_arr)

    class DS(Data.Dataset):
        def __init__(self, p, h, t, l):
            super().__init__()
            self.p, self.h, self.t, self.l = p, h, t, l
        def __len__(self):
            return self.p.shape[0]
        def __getitem__(self, i):
            return self.p[i], self.h[i], self.t[i], self.l[i]

    loader = Data.DataLoader(DS(pep_t, hla_t, tcr_t, lbl_t), batch_size=1024, shuffle=False, num_workers=0)

    # Model params (from PISTE)
    window_map = {'random': '3', 'unipep': '1', 'reftcr': '2'}
    model = Transformer(
        device=device, vocab_size=21, d_model=64, e_layers=3, d=64, n_heads=9,
        sigma=1, window_threshold=window_map.get(model_name, '3'), d_ff=512,
        interact_layers=1, tgt_len=tgt_len, hla_max_len=34, d_layers=1
    ).to(device)

    ckpt = os.path.join(PISTE_CODE, 'checkpoints', model_name, 'exp0.pkl')
    model.load_state_dict(torch.load(ckpt, map_location=device))
    model.eval()

    # Predict
    all_lbls, all_probs = [], []
    with torch.no_grad():
        for bp, bh, bt, bl in loader:
            bp, bh, bt = bp.to(device), bh.to(device), bt.to(device)
            out, _, _ = model(bp, bh, bt)
            probs = F.softmax(out, dim=1)[:, 1].cpu().numpy()
            all_lbls.extend(bl.numpy())
            all_probs.extend(probs)

    y_true, y_prob = np.array(all_lbls), np.array(all_probs)

    # Metrics
    roc_auc = roc_auc_score(y_true, y_prob)
    prec, reca, _ = precision_recall_curve(y_true, y_prob)
    aupr = auc(reca, prec)
    n_pos = int(sum(y_true))
    sorted_idx = sorted(range(len(y_prob)), key=lambda i: y_prob[i], reverse=True)
    sorted_true = [y_true[i] for i in sorted_idx]
    ppvn = sum(sorted_true[:n_pos]) / n_pos if n_pos > 0 else 0

    return {
        'AUC': round(roc_auc, 4), 'AUPR': round(aupr, 4),
        'PPVN': round(ppvn, 4), 'n_pos': n_pos, 'n_neg': len(y_true) - n_pos
    }


def main():
    datasets = {
        'score_gt_0': os.path.join(VDJDB_DIR, 'vdjdb_filtered_score_gt_0_clean.tsv'),
        'score_ge_2': os.path.join(VDJDB_DIR, 'vdjdb_filtered_score_ge_2_clean.tsv'),
    }

    results = {}

    for dname, fpath in datasets.items():
        print(f"\n{'='*70}")
        print(f"Dataset: {dname}")
        print(f"{'='*70}")

        pos_df = load_vdjdb(fpath)
        print(f"Loaded {len(pos_df)} positive triples")

        # Check HLA coverage
        hla_db = set(pd.read_csv(HLA_SEQ_FILE)['HLA_type'])
        unmatched = set(pos_df['HLA_type']) - hla_db
        if unmatched:
            print(f"  Dropping {sum(pos_df['HLA_type'].isin(unmatched))} triples with unmatched HLA: {unmatched}")
            pos_df = pos_df[~pos_df['HLA_type'].isin(unmatched)]
            print(f"  Remaining: {len(pos_df)}")

        if len(pos_df) < 5:
            print(f"  Too few samples, skipping")
            continue

        results[dname] = {}

        for model_name in ['random', 'unipep', 'reftcr']:
            print(f"\n  Model: PISTE-{model_name}")
            try:
                metrics = run_prediction(pos_df, model_name, neg_ratio=10)
                results[dname][model_name] = metrics
                print(f"    AUC={metrics['AUC']:.4f}  AUPR={metrics['AUPR']:.4f}  PPVN={metrics['PPVN']:.4f}")
                print(f"    (n_pos={metrics['n_pos']}, n_neg={metrics['n_neg']})")
            except Exception as e:
                print(f"    ERROR: {e}")
                import traceback
                traceback.print_exc()

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'Dataset':<15} {'Model':<10} {'AUC':<10} {'AUPR':<10} {'PPVN':<10} {'Pos':<8} {'Neg':<8}")
    print("-" * 71)
    for dname in results:
        for mname in results[dname]:
            m = results[dname][mname]
            print(f"{dname:<15} {mname:<10} {m['AUC']:<10.4f} {m['AUPR']:<10.4f} {m['PPVN']:<10.4f} {m['n_pos']:<8} {m['n_neg']:<8}")

    # Save summary
    summary = []
    for dname in results:
        for mname in results[dname]:
            r = results[dname][mname]
            r['dataset'] = dname
            r['model'] = mname
            summary.append(r)
    pd.DataFrame(summary).to_csv(os.path.join(OUTPUT_DIR, 'summary.csv'), index=False)
    print(f"\nSaved to {OUTPUT_DIR}/summary.csv")


if __name__ == '__main__':
    main()
