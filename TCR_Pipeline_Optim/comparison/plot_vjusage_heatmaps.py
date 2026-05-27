"""
plot_vjusage_heatmaps.py — 从 VDJTools PlotFancyVJUsage 的 TXT 输出生成 V-J pairing 热图
"""
import os, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

BASE = "/Users/edwardchan/Desktop/TCR"
VJUSAGE_DIR = os.path.join(BASE, "vjusage")
OUT_DIR = os.path.join(BASE, "comparison", "figures_preset")

SHORT_IDS = ["S044", "S045", "S046", "S047", "S048"]
PRESETS = ["generic-amplicon", "rna-seq"]
PRESET_COLORS = {"generic-amplicon": "Blues", "rna-seq": "Reds"}

os.makedirs(OUT_DIR, exist_ok=True)

# Find all TXT files
all_files = sorted([f for f in os.listdir(VJUSAGE_DIR) if f.endswith(".txt")])
print(f"Found {len(all_files)} VJ usage TXT files")

# Parse one file to understand structure
sample_file = os.path.join(VJUSAGE_DIR, all_files[0])
vj_data = {}

for fname in all_files:
    # Parse sample and preset from filename
    fpath = os.path.join(VJUSAGE_DIR, fname)
    sid = None
    for s in SHORT_IDS:
        if fname.startswith(s):
            sid = s
            break
    preset = "generic-amplicon" if "generic-amplicon" in fname else "rna-seq"

    # Read V-J matrix
    df = pd.read_csv(fpath, sep="\t", index_col=0)
    df.index.name = "J"

    key = f"{sid}_{preset}"
    vj_data[key] = df

print(f"Parsed {len(vj_data)} matrices")

# Get unified V and J gene lists across all samples
all_v = set()
all_j = set()
for df in vj_data.values():
    all_v.update(df.columns)
    all_j.update(df.index)

# Sort genes naturally
def gene_sort_key(g):
    m = re.match(r'TRB([VDJ])(\d+)-?(\d*)', g)
    if m:
        return (m.group(1), int(m.group(2)), int(m.group(3)) if m.group(3) else 0)
    return (g, 0, 0)

v_genes = sorted(all_v, key=gene_sort_key)
j_genes = sorted(all_j, key=gene_sort_key)

print(f"V genes: {len(v_genes)}, J genes: {len(j_genes)}")

# ============================================================
# Figure 1: Per-sample amp vs rna side-by-side (5 samples)
# ============================================================
fig, axes = plt.subplots(5, 2, figsize=(28, 35))

# Find global vmax for consistent color scale (99th percentile to avoid outlier domination)
all_vals = []
for df in vj_data.values():
    all_vals.extend(df.values.flatten())
vmax = np.percentile(all_vals, 99)

for row_idx, sid in enumerate(SHORT_IDS):
    for col_idx, preset in enumerate(PRESETS):
        ax = axes[row_idx, col_idx]
        key = f"{sid}_{preset}"
        if key not in vj_data:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        df = vj_data[key]
        # Reindex to unified gene list, fill missing with 0
        mat = df.reindex(index=j_genes, columns=v_genes).fillna(0)

        cmap = plt.cm.Blues if preset == "generic-amplicon" else plt.cm.Reds

        im = ax.imshow(mat.values, aspect="auto", cmap=cmap,
                       norm=mcolors.LogNorm(vmin=max(vmax*1e-5, mat.values[mat.values > 0].min()),
                                            vmax=max(vmax, 1e-6)),
                       interpolation="none")

        ax.set_xticks(range(len(v_genes)))
        ax.set_xticklabels(v_genes, rotation=90, fontsize=5)
        ax.set_yticks(range(len(j_genes)))
        ax.set_yticklabels(j_genes, fontsize=6)

        title = f"{sid} — {preset}"
        ax.set_title(title, fontsize=10, fontweight="bold")
        plt.colorbar(im, ax=ax, shrink=0.8, label="Frequency")

fig.suptitle("V-J Pairing Heatmaps: generic-amplicon vs rna-seq (log scale, minCount≥2)",
             fontsize=14, fontweight="bold", y=0.998)
plt.tight_layout()
outpath = os.path.join(OUT_DIR, "P14_vj_pairing_all_samples.png")
fig.savefig(outpath, dpi=120, bbox_inches="tight")
print(f"Figure saved: {outpath}")
plt.close()

# ============================================================
# Figure 2: Delta heatmaps (amp - rna) for each sample
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(22, 14))
axes = axes.flatten()

for idx, sid in enumerate(SHORT_IDS):
    ax = axes[idx]
    key_amp = f"{sid}_generic-amplicon"
    key_rna = f"{sid}_rna-seq"

    if key_amp not in vj_data or key_rna not in vj_data:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        continue

    mat_amp = vj_data[key_amp].reindex(index=j_genes, columns=v_genes).fillna(0)
    mat_rna = vj_data[key_rna].reindex(index=j_genes, columns=v_genes).fillna(0)
    mat_delta = mat_amp.values - mat_rna.values

    # Symmetric colormap for delta
    vlim = max(abs(mat_delta.min()), abs(mat_delta.max()))
    # Use 99th percentile to avoid outliers saturating
    flat_deltas = np.abs(mat_delta).flatten()
    vlim_p99 = np.percentile(flat_deltas, 99) if len(flat_deltas) > 0 else vlim

    im = ax.imshow(mat_delta, aspect="auto", cmap="RdBu_r",
                   vmin=-vlim_p99, vmax=vlim_p99, interpolation="none")

    ax.set_xticks(range(len(v_genes)))
    ax.set_xticklabels(v_genes, rotation=90, fontsize=5)
    ax.set_yticks(range(len(j_genes)))
    ax.set_yticklabels(j_genes, fontsize=6)
    ax.set_title(f"{sid}: amp − rna", fontsize=10, fontweight="bold")
    plt.colorbar(im, ax=ax, shrink=0.8, label="Δ Frequency")

axes[-1].axis("off")
fig.suptitle("V-J Pairing Differences: generic-amplicon − rna-seq (Blue=amp higher, Red=rna higher)",
             fontsize=13, fontweight="bold")
plt.tight_layout()
outpath = os.path.join(OUT_DIR, "P14_vj_pairing_delta.png")
fig.savefig(outpath, dpi=120, bbox_inches="tight")
print(f"Figure saved: {outpath}")
plt.close()

# ============================================================
# Figure 3: V-J correlation scatter (amp vs rna, per V-J pair)
# ============================================================
fig, ax = plt.subplots(figsize=(10, 10))

all_amp_vals = []
all_rna_vals = []

for sid in SHORT_IDS:
    key_amp = f"{sid}_generic-amplicon"
    key_rna = f"{sid}_rna-seq"
    if key_amp in vj_data and key_rna in vj_data:
        mat_amp = vj_data[key_amp].reindex(index=j_genes, columns=v_genes).fillna(0)
        mat_rna = vj_data[key_rna].reindex(index=j_genes, columns=v_genes).fillna(0)
        all_amp_vals.extend(mat_amp.values.flatten())
        all_rna_vals.extend(mat_rna.values.flatten())

all_amp_vals = np.array(all_amp_vals)
all_rna_vals = np.array(all_rna_vals)

# Remove zero-zero pairs for clarity
mask = (all_amp_vals > 0) | (all_rna_vals > 0)
amp_filt = all_amp_vals[mask]
rna_filt = all_rna_vals[mask]

# Log transform for visualization (add small offset for zeros)
offset = min(amp_filt[amp_filt > 0].min(), rna_filt[rna_filt > 0].min()) / 10
amp_log = np.log10(amp_filt + offset)
rna_log = np.log10(rna_filt + offset)

ax.scatter(amp_log, rna_log, s=5, alpha=0.5, c="#2196F3", edgecolors="none")

# Correlation
from scipy import stats
r, p = stats.pearsonr(amp_log, rna_log)
ax.plot([amp_log.min(), amp_log.max()], [amp_log.min(), amp_log.max()], 'k--', alpha=0.3)

ax.set_xlabel("generic-amplicon log10(freq)", fontsize=11)
ax.set_ylabel("rna-seq log10(freq)", fontsize=11)
ax.set_title(f"V-J Pairing Frequency Correlation (log scale)\nr = {r:.6f}, p = {p:.2e}",
             fontsize=12, fontweight="bold")
ax.grid(True, alpha=0.3)

outpath = os.path.join(OUT_DIR, "P14_vj_pairing_correlation.png")
fig.savefig(outpath, dpi=120, bbox_inches="tight")
print(f"Figure saved: {outpath}")
plt.close()

print("\nDone. All V-J pairing figures generated.")
