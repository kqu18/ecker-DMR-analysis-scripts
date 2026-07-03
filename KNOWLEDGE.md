# Project Knowledge — Arabidopsis Single-Cell Methylation × RNA

> Self-contained context document. Paste into claude.ai (as a Project knowledge file or first chat message) to discuss this research without launching Claude Code.
>
> Snapshot date: 2026-06-27. Code/path claims are point-in-time; verify before relying on specifics.

---

## 1. Project Goal

Single-cell methylation × RNA expression analysis for *Arabidopsis thaliana* leaf tissue. The goal is to identify **genes** (and **DMWs** — Differentially Methylated Windows) whose **gene-body CG methylation** is **negatively correlated with expression** across **17 cell type clusters**.

Three genotypes are involved: **col** (wild type), **rdd** (RNA-directed DNA methylation mutant), and **met** (T-MET / DNA-methyltransferase mutant).

A manuscript is in preparation (Salk Institute affiliation; see git history for the WLS-correlation methods commits).

---

## 2. Cell Cluster Reference

17 clusters from Arabidopsis leaf, grouped Mesophyll → Epidermis → Other:

```python
cluster_mapping = {
    3: '3 Early M',   0: '0 Mid M',     1: '1 Late M',     7: '7 Expanding M',   10: '10 Senescent M',
    5: '5 Early E',   4: '4 Abaxial E', 6: '6 Adaxial E',  13: '13 Expanding E', 8: '8 Senescent E',
    14: '14 Guard cell', 2: '2 Vasculature', 11: '11 PPP', 9: '9 Phloem',
    16: '16 Myrosinase', 12: '12 S phase', 15: '15 G2M phase'
}
```

Heatmap column order: Mesophyll (0–4) → Epidermis (5–9) → Other (10–16). Vertical dividers mark the Meso/Epi boundary at column 5 and the Epi/Other boundary at column 10. Dotted lines mark senescent clusters within each block.

`cluster_mapping` is sometimes accessed with str keys, sometimes int — keep both: `cluster_mapping.update({str(k): v for k, v in cluster_mapping.items()})`.

---

## 3. Six-Stage Pipeline Arc

Every analysis fits into one of these stages. Understanding the arc tells you which existing notebook to extend.

### Stage 1 — DMR calling (upstream)

Chi-square test on clustered bisulfite data with per-cluster logit offsets, phi overdispersion correction, and BH-FDR. Produces **3,124 significant CG DMWs** in `chunks_CG_minfilt.gff`. Statistical framework developed in `notebooks/2025/presentable.ipynb`. Code in §7 below.

### Stage 2 — MCDS construction

DMW MCDS built externally at `/ceph/MethDev/pbio/andy/JW/section_4/All_mcds_dmw3_{0,1}.mcds`. Contains **all three genotypes** in one file, labeled only by cell-ID prefix (see §6). Concatenation sanity-checked in `section4/combine_mcds.ipynb`.

### Stage 3 — DMW curation

The 3,124 chunks become the canonical row set after filtering: aggregate coverage ≥ 50, max raw fraction ≥ 0.1 in any cluster. Curation lives in `notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb` and `notebooks/2026/260514_dmw_methylation_matrices_by_genotype.ipynb`.

### Stage 4 — DMW × 17-cluster methylation matrices

For each DMW × cluster:
1. Aggregate `mc` and `cov` counts.
2. Fit global per-cluster logit offsets (`deltas`) to remove bulk methylation differences.
3. Compute residual `obs_p − p0` (observed minus cluster-adjusted expected).

Same algorithm appears in all four DMW notebooks — candidate for factoring into a `dmw_utils.compute_logit_residuals()` helper.

### Stage 5 — WLS correlation against RNA

For each gene (or DMW), fit `RNA_z ~ meth_z` weighted by coverage using `statsmodels.WLS`. Extract slope (pseudo-ρ), p-value, BH-FDR q-value.

- Gene-level: `notebooks/2026/260305_corr_heatmap_refactored.ipynb`
- DMW-level: `notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb` — adds `bioframe.overlap()` to map DMW coordinates to gene bodies before WLS

Output: `data/dmw_gene_wls_stats.tsv` (slope, p, q, gene mapping per DMW). WLS code in §8 below.

### Stage 6 — Heatmaps and downstream catalogs

Three figure families:
- **Top-N anticorrelation heatmaps** — methylation (coolwarm) + RNA (custom colormap) side-by-side, Ward-clustered rows on methylation z-scores clipped to ±2.5.
- **Hypomethylation catalogs** — cluster-specific hypo events at gene level (`260325_identify_hypo.ipynb`) and DMW level (`260501_dmw_identify_hypo.ipynb`).
- **Genotype comparison panels** — col / rdd / met side-by-side residual heatmaps anchored to col-0 row set (`260514_dmw_methylation_matrices_by_genotype.ipynb`).

---

## 4. Data Paths

Input data lives **outside** the repo (`/ceph/MethDev/pbio/kay/`). All notebooks use absolute paths.

| Variable | Path |
|---|---|
| Gene MCDS | `/ceph/MethDev/pbio/hazel/gene_mcds/gene` |
| RNA counts | `/ceph/MethDev/pbio/hazel/count_Arab13_with_TEs.csv` |
| Cell cluster metadata | `/ceph/MethDev/pbio/hazel/merged_cluster_assignments.csv` |
| DMW MCDS (two files) | `/ceph/MethDev/pbio/andy/JW/section_4/All_mcds_dmw3_{0,1}.mcds` |
| Gene annotation BED | `/ceph/MethDev/pbio/andy/JW/section_4/genes.fixed.bed` |
| Raw DMR call input | `/ceph/MethDev/pbio/data/annotated_filtered_col.CG_2.fast.tsv` |

When adapting a notebook to a new dataset, **edit the paths block near the top of Cell 2** (look for `# ---- EDIT PATHS HERE ----`).

### Output conventions

| Type | Location |
|---|---|
| Heatmaps (SVG/PNG) | `/ceph/MethDev/pbio/kay/figures/heatmaps/` |
| Per-cluster gene lists (CSV) | `/ceph/MethDev/pbio/kay/data/` |
| Top-gene text lists | `/ceph/MethDev/pbio/kay/genes/` |
| WLS stats table | `/ceph/MethDev/pbio/kay/data/dmw_gene_wls_stats.tsv` |

---

## 5. Canonical Notebooks

The repo has ~55 notebooks. Most are scratchpads. The paper-worthy core is small:

| Notebook | Role |
|---|---|
| `notebooks/2026/260305_corr_heatmap_refactored.ipynb` | **Canonical gene-level WLS pipeline** |
| `notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb` | **Canonical DMW-level WLS pipeline** (DMW→gene overlap, top-300 heatmap, `dmw_gene_wls_stats.tsv`) |
| `notebooks/2026/260514_dmw_methylation_matrices_by_genotype.ipynb` | **Canonical multi-genotype comparison** (col/rdd/met) |
| `notebooks/2026/260325_identify_hypo.ipynb` | Gene-level hypomethylation catalog (cluster-specific) |
| `notebooks/2026/260501_dmw_identify_hypo.ipynb` | DMW-level hypomethylation catalog |
| `notebooks/2025/presentable.ipynb` | DMR-calling statistical framework |

**Useful supplementary:**
- `notebooks/2026/260318_9conditions.ipynb` — 9-cluster-subset slicing → `rna_meth_*` heatmaps
- `notebooks/2026/260331_30_heatmaps.ipynb` — combinatorial heatmap atlas
- `notebooks/2026/heatmap_hazel_data_refined.ipynb` — Hazel's curated gene-list heatmaps
- `notebooks/2026/260401_refine_WLS.ipynb` + `260406_refine_WLS_heatmap.ipynb` — WLS robustness
- `section4/minfilt_gff_analysis.ipynb` — bulk-ALLC vs MCDS validation
- `section4/combine_mcds.ipynb` — MCDS concatenation utility

**Naming convention:** `YYMMDD_description.ipynb` (2026 phase) or `MMDD_description.ipynb` (2025). The highest-dated file for a topic is the most current; `_copy` and `_refactored` suffixes indicate iterative refinements.

---

## 6. Mutant Genotype Encoding

The DMW MCDS contains **all three genotypes** in one file. Genotype lives in the second underscore-field of `cell_id` (format: `240614_mct_<G>_<plate>_<well>`):

| Prefix field | Genotype | Plate range | Cells |
|---|---|---|---|
| `mct_1_*` | **col** (wild type) | 1_1 – 1_4 | 3,059 |
| `mct_2_*` | **rdd** | 2_1 – 2_4 | 1,480 |
| `mct_3_*` | **met** (T-MET) | 3_1 – 3_4 | 1,532 |

Total: 6,071 cells across 12 plates.

To split by genotype:
```python
genotype_map = {"1": "col", "2": "rdd", "3": "met"}
genotypes = pd.Series(mcds["cell"].values).str.split("_").str[2].map(genotype_map)
```

The MCDS has **no genotype coordinate** — only `cell`, `dmw`, `mc_type`, `count_type`. Always derive from the cell-ID prefix.

---

## 7. DMR Calling Code (Stage 1)

Chi-square pipeline producing the DMW set.

```python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import logit, expit
from scipy.stats import chi2
from joblib import Parallel, delayed

df = pd.read_csv('/ceph/MethDev/pbio/data/annotated_filtered_col.CG_2.fast.tsv', sep='\t')
df = df.dropna(subset=['score_masked'])
```

`c` = methylated reads, `t` = unmethylated reads, `score` = methylation ratio, `cluster` = cell cluster ID.

### Global cluster offsets

```python
def fit_offsets(df_all):
    agg = df_all.groupby('cluster')[['c','t']].sum()
    M = agg.sum(axis=1)
    p = agg['c'] / M
    d = logit(np.clip(p, 1e-6, 1-1e-6))
    return d - np.average(d, weights=M)
```

### Window-wise chi-square

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

### BH FDR

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

### Driver

```python
deltas = fit_offsets(df)
stats  = window_stats(df, deltas, tau=20.0)
stats['qval'] = bh_fdr(stats['pval'].to_numpy())
stats.to_csv('/ceph/MethDev/pbio/data/x2_CG_stat_all_std.csv', index=False)
```

`delta_max` = within-window methylation deviation (effect size).

---

## 8. WLS Correlation Code (Stage 5)

Canonical pattern shared by all correlation notebooks:

1. **Load**: `MCDS.open(mcds_paths=..., var_dim='gene'|'dmw')`, RNA CSV (transposed), cluster metadata.
2. **RNA normalization**: raw counts → CPM → log1p → mean per cluster → `cluster_rna_matrix` (genes × clusters).
3. **Methylation residuals** (`filtered_meth_df` / `adj_meth_matrix`):
   - Aggregate mc/cov per cluster.
   - Fit per-cluster logit offsets `deltas`.
   - Compute `obs_p − p0` residuals.
   - Filter: aggregate coverage ≥ 50, max raw fraction ≥ 0.1 in any cluster.
4. **WLS**: for each gene/DMW, fit `RNA_z ~ meth_z` with `statsmodels.WLS` using coverage as weights. Extract slope (pseudo-ρ) and p-value. Apply BH FDR.
5. **Visualize**: side-by-side methylation (coolwarm) + RNA (custom colormap) heatmaps, rows hierarchically clustered by methylation pattern (Ward linkage, z-scores clipped ±2.5).

DMW pipeline adds a `bioframe.overlap()` step to map DMW coordinates → gene bodies before WLS.

---

## 9. Gotchas

- `rna = pd.read_csv(rna_path).T` — the RNA CSV is genes × cells; transpose immediately on load so rows are cells.
- DMW chromosome names in MCDS are bare integers (`1`–`5`); gene BED uses `chr1`-style. Harmonize before `bioframe.overlap()`:
  ```python
  dmw_bed["chrom"] = "chr" + dmw_bed["chrom"].astype(str)
  ```
- `cluster_mapping` must include both int and str keys — cluster indices come in either type depending on the operation.
- `logit`/`expit` require clipping to `(1e-6, 1-1e-6)` before application to avoid ±∞.
- DMW MCDS has no `genotype` coordinate — derive from cell-ID prefix (§6).
- 2025 notebooks contain several pre-existing broken paths (`Arab10_All_mcds.mcds`, `annotated_filtered_col.CG_2.fast.tsv` at the kay/data/ location, `result/highMCClusterMap.pdf`, `MappingSummary.csv.gz`, `gffs_by_clst_genotype`, `master_gffs/`). These were broken before the 2026-05 repo restructuring, not caused by it.

---

## 10. Repository Structure

```
/ceph/MethDev/pbio/kay/
├── notebooks/
│   ├── 2025/     MMDD-prefixed notebooks from the 2025 analysis phase
│   └── 2026/     YYMMDD-prefixed notebooks from the 2026 analysis phase
├── section4/     Shell scripts + notebooks for MCDS generation
├── data/         Input/output data files (data/outliers/)
├── figures/      Output figures (heatmaps/, dmw_heatmaps/, ...)
├── genes/        Gene list text files
├── CLAUDE.md     Project instructions for Claude Code
├── KNOWLEDGE.md  (this file)
└── README.md
```

### Environment

Notebooks use a conda/mamba environment with `ALLCools`, `statsmodels`, `bioframe`, `scipy`, `seaborn`, `xarray`. If imports fail, activate the right env first.

Run a notebook non-interactively:
```bash
jupyter nbconvert --to notebook --execute <notebook>.ipynb --output <notebook>_out.ipynb
```

---

## 11. Shared Logic Worth Factoring Out

These patterns appear copy-pasted across 4+ notebooks — candidates for a shared utilities module:

- Per-cluster mc/cov aggregation + logit-residual computation
- `bioframe.overlap()` DMW→gene mapping with `chr` prefix harmonization
- `cluster_mapping` dict + `cluster_order` list (Meso → Epi → Other)
- Heatmap helper with Meso/Epi/Other dividers and senescent dotted lines

---

## 12. Open Questions / Future Directions

These are the kinds of things worth discussing in a chat session:

- Consolidating canonical notebooks (§5) into a single master pipeline notebook.
- Extending the col-0 DMW residual pipeline to produce `rdd_dmw_*` and `met_dmw_*` matrices (subset cells before per-cluster aggregation; reuse the 17-cluster IDs from `merged_cluster_assignments.csv`).
- Factoring shared logic (§11) into `dmw_utils.py`.
- Diagnosing the 2025 broken paths (§9) if any of those notebooks need to be revived.
