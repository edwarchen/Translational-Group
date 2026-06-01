#!/usr/bin/env python3
"""
Evaluate PISTE on external VDJdb test sets.

1. Convert VDJdb filtered data to PISTE input format
2. Generate negative samples (random shuffle & unique peptide)
3. Run PISTE predictions
4. Calculate AUC, AUPR, PPVN, etc.
"""

import os
import sys
import numpy as np
import pandas as pd
import random
from collections import Counter

# Paths
PISTE_CODE = '/home/chenya/TCR_MODEL_TEST/PISTE/code'
VDJDB_DIR = '/home/chenya/TCR_MODEL_TEST/data'
HLA_SEQ_FILE = '/home/chenya/TCR_MODEL_TEST/PISTE/data/raw_data/common_hla_sequence.csv'
OUTPUT_DIR = '/home/chenya/TCR_MODEL_TEST/data/piste_eval'

os.makedirs(OUTPUT_DIR, exist_ok=True)

random.seed(42)
np.random.seed(42)


def load_and_prepare_vdjdb(vdjdb_file):
    """
    Load filtered VDJdb data and convert to PISTE format.
    Returns DataFrame with columns: CDR3, MT_pep, HLA_type
    """
    df = pd.read_csv(vdjdb_file, sep='\t')

    # Rename columns to PISTE format
    df = df.rename(columns={
        'cdr3.beta': 'CDR3',
        'peptide': 'MT_pep',
    })

    # Convert HLA to PISTE format: HLA-A*02:01 -> HLA-A02:01
    df['HLA_type'] = df['hla'].str.replace('*', '', regex=False)

    return df[['CDR3', 'MT_pep', 'HLA_type']]


def map_hla_to_db(hla_types, hla_db):
    """
    Map our HLA types to PISTE's HLA DB entries.
    If exact match not found, try mapping 2-field to most common 4-field subtype.
    """
    db_set = set(hla_db['HLA_type'].values)
    mapped = {}
    unmapped = set()

    for hla in set(hla_types):
        if hla in db_set:
            mapped[hla] = hla
        else:
            # Try to find best match: HLA-A02 -> HLA-A02:01
            # Find all DB entries starting with the same prefix
            matching = [x for x in db_set if x.startswith(hla + ':')]
            if matching:
                # Use the first (typically the most common) one
                mapped[hla] = matching[0]
            else:
                unmapped.add(hla)

    return mapped, unmapped


def generate_negatives_random(pos_df, neg_ratio=10):
    """
    Generate negative samples by randomly pairing TCRs with peptides and HLAs.
    """
    negatives = []
    cdrs = pos_df['CDR3'].values
    peps = pos_df['MT_pep'].values
    hlas = pos_df['HLA_type'].values

    n_neg = len(pos_df) * neg_ratio
    pos_triples = set(zip(cdrs, peps, hlas))

    while len(negatives) < n_neg:
        i = random.randint(0, len(cdrs) - 1)
        j = random.randint(0, len(cdrs) - 1)
        k = random.randint(0, len(cdrs) - 1)

        triple = (cdrs[i], peps[j], hlas[k])
        if triple not in pos_triples:
            negatives.append(triple)
            pos_triples.add(triple)  # Avoid duplicate negatives

    neg_df = pd.DataFrame(negatives, columns=['CDR3', 'MT_pep', 'HLA_type'])
    return neg_df


def generate_negatives_unipep(pos_df, neg_ratio=10):
    """
    Generate negative samples with unique peptide sampling.
    Each positive peptide gets paired with random TCRs not associated with it.
    """
    negatives = []
    pos_triples = set(zip(pos_df['CDR3'], pos_df['MT_pep'], pos_df['HLA_type']))

    # Group TCRs by HLA
    hla_to_tcrs = pos_df.groupby('HLA_type')['CDR3'].apply(list).to_dict()
    hla_to_peps = pos_df.groupby('HLA_type')['MT_pep'].apply(list).to_dict()

    for _, row in pos_df.iterrows():
        hla = row['HLA_type']
        pep = row['MT_pep']
        tcr = row['CDR3']

        # Find TCRs not associated with this peptide for this HLA
        other_tcrs = [t for t in hla_to_tcrs[hla] if (t, pep, hla) not in pos_triples][:neg_ratio]

        # If not enough, sample from broader pool
        if len(other_tcrs) < neg_ratio:
            extra_tcrs = [t for t in hla_to_tcrs[hla] if t != tcr and (t, pep, hla) not in pos_triples]
            other_tcrs = (other_tcrs + extra_tcrs)[:neg_ratio]

        for neg_tcr in other_tcrs:
            negatives.append((neg_tcr, pep, hla))

    neg_df = pd.DataFrame(negatives, columns=['CDR3', 'MT_pep', 'HLA_type'])
    neg_df = neg_df.drop_duplicates()

    # Trim to target ratio
    target_n = len(pos_df) * neg_ratio
    if len(neg_df) > target_n:
        neg_df = neg_df.sample(n=target_n, random_state=42)

    return neg_df


def prepare_eval_files(pos_df, eval_name, neg_method='random', neg_ratio=10):
    """
    Prepare evaluation CSV files with positive and negative samples.
    """
    # Generate negatives
    if neg_method == 'random':
        neg_df = generate_negatives_random(pos_df, neg_ratio)
    elif neg_method == 'unipep':
        neg_df = generate_negatives_unipep(pos_df, neg_ratio)
    else:
        raise ValueError(f"Unknown neg_method: {neg_method}")

    # Label
    pos_df_labeled = pos_df.copy()
    pos_df_labeled['Label'] = 1
    neg_df_labeled = neg_df.copy()
    neg_df_labeled['Label'] = 0

    # Combine
    combined = pd.concat([pos_df_labeled, neg_df_labeled], ignore_index=True)
    combined = combined.sample(frac=1, random_state=42).reset_index(drop=True)

    # Save
    out_path = os.path.join(OUTPUT_DIR, f'{eval_name}_{neg_method}.csv')
    combined.to_csv(out_path, index=False)

    print(f"  {eval_name} ({neg_method}): {len(pos_df)} pos + {len(neg_df)} neg -> {out_path}")
    return out_path


def make_data_for_prediction(df, pep_max_len=11, hla_max_len=34, tcr_max_len=30):
    """
    Convert DataFrame to model inputs (same as PISTE's make_data).
    """
    vocab = {'C': 1, 'W': 2, 'V': 3, 'A': 4, 'H': 5, 'T': 6, 'E': 7, 'K': 8,
             'N': 9, 'P': 10, 'I': 11, 'L': 12, 'S': 13, 'D': 14, 'G': 15,
             'Q': 16, 'R': 17, 'Y': 18, 'F': 19, 'M': 20, '-': 0}

    pep_inputs, hla_inputs, tcr_inputs, labels = [], [], [], []

    for _, row in df.iterrows():
        pep = str(row['MT_pep']).ljust(pep_max_len, '-')
        hla = str(row['HLA_sequence']).ljust(hla_max_len, '-')
        tcr = str(row['CDR3']).ljust(tcr_max_len, '-')

        pep_inputs.append([vocab.get(aa, 0) for aa in pep])
        hla_inputs.append([vocab.get(aa, 0) for aa in hla])
        tcr_inputs.append([vocab.get(aa, 0) for aa in tcr])
        labels.append(row['Label'])

    pep_tensor = np.array(pep_inputs, dtype=np.int64)
    hla_tensor = np.array(hla_inputs, dtype=np.int64)
    tcr_tensor = np.array(tcr_inputs, dtype=np.int64)
    label_tensor = np.array(labels)

    return pep_tensor, hla_tensor, tcr_tensor, label_tensor


def run_piste_prediction(eval_csv, model_name='random', gpu=True):
    """
    Run PISTE model on the evaluation CSV.
    Returns DataFrame with predictions.
    """
    import torch
    import torch.nn as nn
    import torch.utils.data as Data
    import torch.nn.functional as F

    import importlib.util
    spec = importlib.util.spec_from_file_location("PISTE",
        os.path.join(PISTE_CODE, "Model", "PISTE.py"))
    piste_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(piste_module)
    Transformer = piste_module.Transformer

    # Parameters (same as PISTE)
    pep_max_len = 11
    hla_max_len = 34
    tcr_max_len = 30
    tgt_len = pep_max_len + hla_max_len + tcr_max_len
    batch_size = 1024
    d_model = 64
    dim = 64
    d_ff = 512
    e_layers = 3
    n_heads = 9
    sigma = 1
    d_layers = 1
    interact_layers = 1
    window_size = '3'
    vocab_size = 21

    device = torch.device("cuda:0" if (gpu and torch.cuda.is_available()) else "cpu")
    print(f"  Using device: {device}")

    # Load data with HLA sequences
    df = pd.read_csv(eval_csv)
    hla_seq = pd.read_csv(HLA_SEQ_FILE)
    df = pd.merge(df, hla_seq, on='HLA_type', how='inner')
    print(f"  After HLA merge: {len(df)} triples (dropped {len(pd.read_csv(eval_csv)) - len(df)})")

    # Make data
    pep_inputs, hla_inputs, tcr_inputs, labels = make_data_for_prediction(
        df, pep_max_len, hla_max_len, tcr_max_len)

    class MyDataSet(Data.Dataset):
        def __init__(self, p, h, t, l):
            super().__init__()
            self.p = p
            self.h = h
            self.t = t
            self.l = l

        def __len__(self):
            return self.p.shape[0]

        def __getitem__(self, idx):
            return self.p[idx], self.h[idx], self.t[idx], self.l[idx]

    dataset = MyDataSet(
        torch.LongTensor(pep_inputs),
        torch.LongTensor(hla_inputs),
        torch.LongTensor(tcr_inputs),
        torch.LongTensor(labels)
    )
    loader = Data.DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    # Load model
    model_dir = os.path.join(PISTE_CODE, 'checkpoints', model_name)
    checkpoint_path = os.path.join(model_dir, 'exp0.pkl')

    class WeightedFocalLoss(nn.Module):
        def __init__(self, alpha=.25, gamma=2):
            super().__init__()
            self.alpha = torch.tensor([alpha, 1 - alpha]).to(device)
            self.gamma = gamma

        def forward(self, inputs, targets):
            BCE_loss = F.cross_entropy(inputs, targets, reduction='none')
            targets = targets.type(torch.long)
            at = self.alpha.gather(0, targets.data.view(-1))
            pt = torch.exp(-BCE_loss)
            F_loss = at * (1 - pt) ** self.gamma * BCE_loss
            return F_loss.mean()

    model = Transformer(
        device=device, vocab_size=vocab_size, d_model=d_model,
        e_layers=e_layers, d=dim, n_heads=n_heads, sigma=sigma,
        window_threshold=window_size, d_ff=d_ff,
        interact_layers=interact_layers, tgt_len=tgt_len,
        hla_max_len=hla_max_len, d_layers=d_layers
    ).to(device)

    model.load_state_dict(torch.load(checkpoint_path, map_location=device))
    model.eval()

    # Predict
    all_labels = []
    all_probs = []

    with torch.no_grad():
        for batch_pep, batch_hla, batch_tcr, batch_labels in loader:
            batch_pep = batch_pep.to(device)
            batch_hla = batch_hla.to(device)
            batch_tcr = batch_tcr.to(device)

            outputs, _, _ = model(batch_pep, batch_hla, batch_tcr)
            probs = nn.Softmax(dim=1)(outputs)[:, 1].cpu().numpy()

            all_labels.extend(batch_labels.numpy())
            all_probs.extend(probs)

    df['prob'] = all_probs
    df['label'] = all_labels

    return df


def calculate_metrics(y_true, y_prob):
    """Calculate AUC, AUPR, PPVN."""
    from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

    roc_auc = roc_auc_score(y_true, y_prob)
    prec, reca, _ = precision_recall_curve(y_true, y_prob)
    aupr = auc(reca, prec)

    # PPVN: Precision among top-n predictions (n = number of positives)
    n_pos = sum(y_true)
    if n_pos > 0:
        sorted_indices = sorted(range(len(y_prob)), key=lambda i: y_prob[i], reverse=True)
        sorted_true = [y_true[i] for i in sorted_indices]
        ppvn = sum(sorted_true[:n_pos]) / n_pos
    else:
        ppvn = 0

    return {
        'AUC': roc_auc,
        'AUPR': aupr,
        'PPVN': ppvn,
        'n_pos': int(n_pos),
        'n_neg': int(len(y_true) - n_pos),
    }


def main():
    # Load VDJdb data
    datasets = {
        'score_gt_0': os.path.join(VDJDB_DIR, 'vdjdb_filtered_score_gt_0_clean.tsv'),
        'score_ge_2': os.path.join(VDJDB_DIR, 'vdjdb_filtered_score_ge_2_clean.tsv'),
    }

    for name, fpath in datasets.items():
        print(f"\n{'='*70}")
        print(f"Preparing data: {name}")
        print(f"{'='*70}")

        pos_df = load_and_prepare_vdjdb(fpath)

        # Map HLA and filter unmapped
        hla_db = pd.read_csv(HLA_SEQ_FILE)
        mapped, unmapped = map_hla_to_db(pos_df['HLA_type'].unique(), hla_db)
        if unmapped:
            print(f"  Warning: {len(unmapped)} HLA types not found in PISTE DB: {unmapped}")
            print(f"    These triples will be dropped.")

        pos_df['HLA_type'] = pos_df['HLA_type'].map(mapped)
        pos_df = pos_df.dropna(subset=['HLA_type'])
        print(f"  After HLA mapping: {len(pos_df)} triples")

        if len(pos_df) < 10:
            print(f"  Too few triples, skipping...")
            continue

        # Prepare eval files for both sampling methods
        for neg_method in ['random']:
            print(f"\n  --- Negative sampling: {neg_method} ---")

            # Create eval CSV
            eval_csv = prepare_eval_files(pos_df, name, neg_method=neg_method, neg_ratio=10)

            # Run all three PISTE models
            for model_name in ['random', 'unipep', 'reftcr']:
                print(f"\n  Running PISTE model: {model_name}")
                try:
                    results_df = run_piste_prediction(eval_csv, model_name=model_name, gpu=True)

                    # Calculate metrics
                    y_true = results_df['label'].values
                    y_prob = results_df['prob'].values
                    metrics = calculate_metrics(y_true, y_prob)

                    print(f"    Results: {metrics}")

                    # Save predictions
                    pred_out = os.path.join(OUTPUT_DIR, f'{name}_{neg_method}_{model_name}_predictions.csv')
                    results_df.to_csv(pred_out, index=False)

                    # Save metrics
                    metrics_out = os.path.join(OUTPUT_DIR, f'{name}_{neg_method}_{model_name}_metrics.txt')
                    with open(metrics_out, 'w') as f:
                        for k, v in metrics.items():
                            f.write(f"{k}: {v}\n")

                    print(f"    Saved to: {pred_out}")

                except Exception as e:
                    print(f"    Error running {model_name}: {e}")
                    import traceback
                    traceback.print_exc()

    print(f"\n{'='*70}")
    print("Evaluation complete!")
    print(f"Results saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
