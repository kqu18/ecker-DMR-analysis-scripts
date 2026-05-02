# Ecker Lab — Gene-Body Methylation × Expression Analysis

Downstream analysis repository for the single-nucleus multi-omic atlas of *Arabidopsis thaliana* rosette leaves (Walker et al., Salk Institute). Identifies genes and differentially methylated windows (DMWs) whose gene-body CG methylation is negatively correlated with RNA expression across 17 cell type clusters.

---

## Scientific Background

This repository supports analysis from a ploidy-resolved single-nucleus atlas of *Arabidopsis* rosette leaves, generated using **snmCT-seq** (single-nucleus methylcytosine and transcriptome sequencing). Nuclei from 21-day long-day rosettes were flow-sorted into ploidy-enriched gates (2C–16C) and profiled for parallel DNA methylation and transcriptome, yielding an atlas of 17 annotated cell types spanning mesophyll, epidermis, vasculature, guard cells, phloem, and cell-cycle states.

**Key findings from the atlas:**

- **Methylation decreases with ploidy** — Global CG, CHG, and CHH methylation levels decline progressively with increasing ploidy in a cell-type- and chromatin-dependent manner.
- **Transposon activation accompanies methylation loss** — Ploidy-associated and cell-type-specific methylation erosion is linked to increased expression of transposon superfamilies (LTR/Gypsy, LTR/Copia, DNA/MuDR, LINE/L1, etc.).
- **Phloem companion cells fail to remethylate** — 4C phloem nuclei show near-complete loss of CG methylation compared to 2C counterparts, a failure not seen in other cell types. Single-molecule profiles confirm substantial CG erosion after replication.
- **Small RNA pathway activity in phloem** — RdDM (RNA-directed DNA methylation) small RNA biogenesis factors are enriched in phloem companion cells, whereas core DNA methylation maintenance machinery shows lower expression, suggesting an alternative methylation maintenance strategy.

**This repository's question:** Within this atlas, which genes show cell-type-specific gene-body CG methylation that is negatively correlated with their RNA expression? Identifying these genes connects the observed epigenomic remodelling to transcriptional regulation across the 17 cluster atlas.

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
