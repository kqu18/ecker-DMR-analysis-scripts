# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell methylation × RNA expression analysis for *Arabidopsis thaliana* leaf tissue. The goal is to identify genes (and DMWs — Differentially Methylated Windows) whose gene-body CG methylation is negatively correlated with expression across 17 cell type clusters.

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

Input data lives **outside** this directory. All notebooks reference it via relative paths:

| Variable | Path |
|---|---|
| Gene MCDS | `../../hazel/gene_mcds/gene` |
| RNA counts | `../../hazel/count_Arab13_with_TEs.csv` |
| Cell cluster metadata | `../../hazel/merged_cluster_assignments.csv` |
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

The DMW notebook (`260428_dmw_corr_heatmap_claude.ipynb`) adds a gene-overlap step using `bioframe.overlap()` to map DMW coordinates to gene bodies before running the same WLS.

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
| Heatmaps (SVG/PNG) | `./figures/heatmaps/` |
| Per-cluster gene lists (CSV) | `./data/` |
| Top-gene text lists | `./genes/` |
| WLS stats table | `dmw_gene_wls_stats.tsv` (root) |

## Notebook Naming

Notebooks are named `YYMMDD_description.ipynb`. The most current version of each analysis is the highest-dated file; `_copy` and `_refactored` suffixes indicate iterative refinements within a date.

Key notebooks:
- `260305_corr_heatmap_refactored.ipynb` — canonical gene-level WLS pipeline
- `260428_dmw_corr_heatmap_claude.ipynb` — DMW-level pipeline (well-documented, good template for new analyses)
- `260325_identify_hypo.ipynb` — identifies cluster-specific hypomethylated genes using residual + z-score thresholds

## Common Gotchas

- `rna = pd.read_csv(rna_path).T` — the RNA CSV is genes × cells; transpose immediately on load so rows are cells
- DMW chromosome names in MCDS are bare integers (1–5); gene BED uses `chr1`-style. Harmonize with `dmw_bed["chrom"] = "chr" + dmw_bed["chrom"].astype(str)` before `bioframe.overlap()`
- `cluster_mapping` must include both integer and string keys (`dict.update({str(k): v ...})`) because cluster indices can be either type depending on the operation
- `logit`/`expit` require clipping to `(1e-6, 1-1e-6)` before application to avoid −∞/+∞
