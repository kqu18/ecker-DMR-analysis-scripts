# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell methylation × RNA expression analysis for *Arabidopsis thaliana* leaf tissue. The goal is to identify genes (and DMWs — Differentially Methylated Windows) whose gene-body CG methylation is negatively correlated with expression across 17 cell type clusters.

## Repository Structure

```
.
├── notebooks/
│   ├── 2025/     MMDD-prefixed notebooks from the 2025 analysis phase
│   └── 2026/     YYMMDD-prefixed notebooks from the 2026 analysis phase
├── section4/     Shell scripts and notebooks for MCDS generation
├── data/         All input/output data files and outlier CSVs (data/outliers/)
├── figures/      All output figures (heatmaps/, dmw_heatmaps/, etc.)
├── genes/        Gene list text files
├── CLAUDE.md     (this file)
└── README.md
```

## Running Notebooks

All work is done in Jupyter notebooks. Launch with:
```bash
jupyter notebook
# or
jupyter lab
```

Run a single notebook non-interactively:
```bash
jupyter nbconvert --to notebook --execute <notebook>.ipynb --output <notebook>_out.ipynb
```

The kernel uses a conda/mamba environment that includes `ALLCools`, `statsmodels`, `bioframe`, `scipy`, `seaborn`, and `xarray`. If a notebook fails to import `ALLCools`, activate the correct environment first.

## Data Paths

Input data lives **outside** this directory. All notebooks reference it via absolute paths:

| Variable | Path |
|---|---|
| Gene MCDS | `/ceph/MethDev/pbio/hazel/gene_mcds/gene` |
| RNA counts | `/ceph/MethDev/pbio/hazel/count_Arab13_with_TEs.csv` |
| Cell cluster metadata | `/ceph/MethDev/pbio/hazel/merged_cluster_assignments.csv` |
| DMW MCDS (two files) | `/ceph/MethDev/pbio/andy/JW/section_4/All_mcds_dmw3_{0,1}.mcds` |
| Gene annotation BED | `/ceph/MethDev/pbio/andy/JW/section_4/genes.fixed.bed` |

When adapting a notebook to a new dataset, **edit the paths block near the top of Cell 2** (look for `# ---- EDIT PATHS HERE ----`).

## Core Analysis Pipeline

Every notebook follows this pattern:

1. **Load**: `MCDS.open(mcds_paths=..., var_dim='gene'|'dmw')` + RNA CSV (transposed on load) + cluster metadata CSV
2. **RNA normalization**: raw counts → CPM → log1p → mean per cluster → `cluster_rna_matrix` (genes × clusters)
3. **Methylation residuals** (`filtered_meth_df` / `adj_meth_matrix`):
   - Aggregate `mc` and `cov` counts per cluster
   - Fit global per-cluster logit offsets (`deltas`) to remove bulk methylation differences
   - Compute `obs_p − p0` residuals (observed minus cluster-adjusted expected)
   - Filter: coverage ≥ 50 (aggregate), max raw fraction ≥ 0.1 in at least one cluster
4. **WLS correlation**: for each gene/DMW, fit `RNA_z ~ meth_z` using `statsmodels.WLS` with coverage as weights; extract slope (pseudo-ρ) and p-value; apply BH FDR correction
5. **Visualize**: side-by-side heatmaps — methylation (coolwarm) and RNA (custom colormap) — with rows hierarchically clustered by methylation pattern (Ward linkage on methylation z-scores clipped to ±2.5)

The DMW notebook (`notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb`) adds a gene-overlap step using `bioframe.overlap()` to map DMW coordinates to gene bodies before running the same WLS.

## Cell Cluster Reference

17 clusters from Arabidopsis leaf. Numeric IDs map to:

```python
cluster_mapping = {
    3: '3 Early M', 0: '0 Mid M', 1: '1 Late M', 7: '7 Expanding M', 10: '10 Senescent M',
    5: '5 Early E', 4: '4 Abaxial E', 6: '6 Adaxial E', 13: '13 Expanding E', 8: '8 Senescent E',
    14: '14 Guard cell', 2: '2 Vasculature', 11: '11 PPP', 9: '9 Phloem',
    16: '16 Myrosinase', 12: '12 S phase', 15: '15 G2M phase'
}
```

Heatmap column order: Mesophyll (0–4) → Epidermis (5–9) → Other cell types (10–16). Vertical dividers mark the Meso/Epi boundary at column 5 and the Epi/Other boundary at column 10.

## Output Conventions

| Type | Location |
|---|---|
| Heatmaps (SVG/PNG) | `/ceph/MethDev/pbio/kay/figures/heatmaps/` |
| Per-cluster gene lists (CSV) | `/ceph/MethDev/pbio/kay/data/` |
| Top-gene text lists | `/ceph/MethDev/pbio/kay/genes/` |
| WLS stats table | `/ceph/MethDev/pbio/kay/data/dmw_gene_wls_stats.tsv` |

## Notebook Naming

Notebooks are named `YYMMDD_description.ipynb` (2026 phase) or `MMDD_description.ipynb` (2025 phase). The most current version of each analysis is the highest-dated file; `_copy` and `_refactored` suffixes indicate iterative refinements within a date.

Key notebooks:
- `notebooks/2026/260305_corr_heatmap_refactored.ipynb` — canonical gene-level WLS pipeline
- `notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb` — DMW-level pipeline (well-documented, good template for new analyses)
- `notebooks/2026/260325_identify_hypo.ipynb` — identifies cluster-specific hypomethylated genes using residual + z-score thresholds

## Common Gotchas

- `rna = pd.read_csv(rna_path).T` — the RNA CSV is genes × cells; transpose immediately on load so rows are cells
- DMW chromosome names in MCDS are bare integers (1–5); gene BED uses `chr1`-style. Harmonize with `dmw_bed["chrom"] = "chr" + dmw_bed["chrom"].astype(str)` before `bioframe.overlap()`
- `cluster_mapping` must include both integer and string keys (`dict.update({str(k): v ...})`) because cluster indices can be either type depending on the operation
- `logit`/`expit` require clipping to `(1e-6, 1-1e-6)` before application to avoid −∞/+∞

---

## DMR Calling Methodology

The following documents the chi-square pipeline used to call Differentially Methylated Regions (DMRs) from clustered bisulfite sequencing data.

### Setup: Import Libraries

```python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import logit, expit
from scipy.stats import chi2
from joblib import Parallel, delayed
```

### Data Loading

Data is loaded from a tab-separated file. Each row represents a genomic location within a cell cluster.

- **`c`**: methylated reads; **`t`**: unmethylated reads; **`score`**: methylation ratio; **`cluster`**: cell cluster ID.

```python
df = pd.read_csv('/ceph/MethDev/pbio/data/annotated_filtered_col.CG_2.fast.tsv', sep='\t')
df = df.dropna(subset=['score_masked'])
```

### Calculate Global Cluster Offsets (`fit_offsets`)

Calculates a baseline methylation level per cluster across the whole genome. Used to establish the null hypothesis.

```python
def fit_offsets(df_all):
    agg = df_all.groupby('cluster')[['c','t']].sum()
    M = agg.sum(axis=1)
    p = agg['c'] / M
    d = logit(np.clip(p, 1e-6, 1-1e-6))
    return d - np.average(d, weights=M)
```

### Window-wise Chi-Square Test (`window_stats`)

For each genomic window, tests whether observed methylation deviates significantly from the cluster-adjusted null.

```python
def window_stats(df_all, deltas, tau=20.0, eps=1e-6):
    out = []
    for (chr_, start, end), g in df_all.groupby(['chr','start','end']):
        c = g['c'].to_numpy(); t = g['t'].to_numpy(); m = c + t
        k = g['cluster'].to_numpy()
        keep = m >= 5
        if keep.sum() < 2:
            continue
        c, m, k = c[keep], m[keep], k[keep]
        pbar = c.sum() / m.sum()
        p0 = expit(logit(np.clip(pbar, 1e-6, 1-1e-6)) + deltas.loc[k].to_numpy())
        E = m * p0
        Var = m * p0 * (1 - p0) + eps
        X2 = ((c - E)**2 / Var).sum()
        dfree = len(c) - 1
        p_hat = c / m
        dev = np.abs(p_hat - p0)
        dmax = float(dev.max() - dev.min())
        hi = int(k[p0.argmax()])
        lo = int(k[p0.argmin()])
        out.append((chr_, start, end, X2, m.sum(), dfree, dmax, hi, lo))
    res = pd.DataFrame(out, columns=['chr','start','end','X2','∑m','df','delta_max','hi_cluster','lo_cluster'])
    phi = np.median(res['X2'] / np.maximum(res['df'], 1)) if len(res) else 1.0
    res['pval'] = 1 - chi2.cdf(res['X2'] / max(phi, 1e-6), res['df'])
    res['phi'] = phi
    return res
```

### BH FDR Correction

```python
def bh_fdr(p):
    if len(p) == 0:
        return p
    r = np.argsort(p)
    ranks = np.empty_like(r); ranks[r] = np.arange(1, len(p) + 1)
    q = p * len(p) / np.maximum(ranks, 1)
    q_sorted = np.minimum.accumulate(np.sort(q)[::-1])[::-1]
    out = np.empty_like(q_sorted); out[r] = q_sorted
    return np.clip(out, 0, 1)
```

### Running the Pipeline

```python
deltas = fit_offsets(df)
stats = window_stats(df, deltas, tau=20.0)
stats['qval'] = bh_fdr(stats['pval'].to_numpy())
stats.to_csv('/ceph/MethDev/pbio/data/x2_CG_stat_all_std.csv', index=False)
```

`delta_max` in the output is the maximum within-window methylation deviation — a measure of effect size.
