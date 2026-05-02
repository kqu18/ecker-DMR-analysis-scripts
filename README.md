# Ecker DMR Analysis Scripts

Analysis notebooks and scripts for single-cell methylation × RNA expression analysis in *Arabidopsis thaliana* leaf tissue. Identifies genes and differentially methylated windows (DMWs) whose gene-body CG methylation is negatively correlated with expression across 17 cell type clusters.

---

## Repository Structure

```
.
├── notebooks/
│   ├── 2025/     Early-phase analysis notebooks (DMG calling, EDA, correlation)
│   └── 2026/     Current analysis notebooks (WLS correlation, heatmaps, DMW pipeline)
├── section4/     Shell scripts and notebooks for MCDS generation (allcools)
├── data/         Input/output data files and outlier CSVs (data/outliers/)
├── figures/      Output figures (heatmaps/, dmw_heatmaps/, etc.)
├── genes/        Gene list text files
├── CLAUDE.md     Project guide and DMR methodology reference
└── README.md     (this file)
```

---

## Key Notebooks

| Notebook | Description |
|---|---|
| `notebooks/2026/260305_corr_heatmap_refactored.ipynb` | Canonical gene-level WLS correlation pipeline |
| `notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb` | DMW-level pipeline; template for new analyses |
| `notebooks/2026/260325_identify_hypo.ipynb` | Identifies cluster-specific hypomethylated genes |
| `notebooks/2026/260501_dmw_identify_hypo.ipynb` | DMW-level hypomethylation analysis |

---

## Dependencies

The analysis environment requires:

- [ALLCools](https://github.com/lhqing/ALLCools) — MCDS data format and methylation utilities
- `statsmodels`, `scipy`, `pandas`, `numpy`
- `bioframe` — genomic interval operations
- `seaborn`, `matplotlib`

Set up with conda/mamba:

```bash
conda activate allcools_kay
jupyter lab
```

---

## Data

Input data lives outside this repository on the cluster:

| Data | Path |
|---|---|
| Gene MCDS | `/ceph/MethDev/pbio/hazel/gene_mcds/gene` |
| RNA counts | `/ceph/MethDev/pbio/hazel/count_Arab13_with_TEs.csv` |
| Cell cluster metadata | `/ceph/MethDev/pbio/hazel/merged_cluster_assignments.csv` |
| DMW MCDS | `/ceph/MethDev/pbio/andy/JW/section_4/All_mcds_dmw3_{0,1}.mcds` |

---

## Contact

For questions, open a GitHub Issue or contact the corresponding author.
