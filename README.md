# Ecker Lab — Gene-Body Methylation × Expression Analysis

Downstream analysis repository for the single-nucleus multi-omic atlas of *Arabidopsis thaliana* rosette leaves (Walker et al., Salk Institute). Identifies genes and differentially methylated windows (DMWs) whose gene-body CG methylation is negatively correlated with RNA expression across 17 cell type clusters.

---

## Scientific Background

### Concepts for the data-oriented reader

| Term | Plain-language meaning |
|---|---|
| **DNA methylation** | A chemical tag (a methyl group) added to cytosine bases in DNA. It does not change the DNA sequence itself, but it can silence genes or alter how actively they are read. Think of it as a sticky note on a specific word in a book that tells the cell "read this less." |
| **CG / CHG / CHH methylation** | Three sequence contexts where the methyl tag can occur. CG (cytosine followed by guanine) is the most studied and most stable context; it is also the focus of gene-body methylation research. |
| **Gene-body methylation** | Methylation found *within* a gene's coding region (not at its promoter). In plants, higher gene-body CG methylation tends to correlate with *lower* expression — the signal this repo specifically looks for. |
| **Ploidy / endoreduplication** | Normal plant cells are diploid (2 copies of each chromosome, called 2C). Plants often allow cells to duplicate their DNA without dividing, producing 4C, 8C, or 16C cells. This is called endoreduplication. Many leaf cell types do this as part of normal growth. |
| **snmCT-seq** | A sequencing technique that simultaneously reads DNA methylation and gene expression (RNA) from the same single nucleus — think of it as running ATAC-seq and RNA-seq on the same cell, but for methylation instead of chromatin accessibility. |
| **Cell clusters / cell types** | Cells grouped by similarity in their gene expression profiles (same idea as scRNA-seq clustering). Each cluster corresponds to a distinct cell type or developmental state. This atlas identifies 17 such clusters. |
| **Transposons** | Segments of DNA that can copy and insert themselves elsewhere in the genome ("jumping genes"). Normally silenced by methylation; when methylation is lost, they can become active, potentially destabilizing the genome. |
| **RdDM (RNA-directed DNA methylation)** | A pathway that uses small RNAs to guide methylation machinery back to specific genomic locations — a mechanism for re-establishing methylation after it has been lost. |

### The study

This repository supports analysis from a ploidy-resolved single-nucleus atlas of *Arabidopsis* rosette leaves (Walker et al., Salk Institute), generated using **snmCT-seq**. Nuclei from 21-day-old rosettes were flow-sorted by ploidy (2C–16C) and profiled for both DNA methylation and RNA expression, producing an atlas of 17 annotated cell types spanning mesophyll, epidermis, vasculature, guard cells, phloem, and cell-cycle states.

**Key findings from the atlas:**

- **Methylation decreases with ploidy** — As cells undergo more rounds of endoreduplication (2C → 4C → 8C → 16C), their CG, CHG, and CHH methylation levels progressively decline. This varies by cell type and by which parts of the genome are measured.
- **Transposon activation accompanies methylation loss** — The methylation erosion is not just a passive bookkeeping change: transposons that were kept silent by methylation become active in high-ploidy cells, a potential threat to genome integrity.
- **Phloem companion cells fail to remethylate** — Most cell types partially recover methylation after each replication round. Phloem companion cells are an exception: their 4C nuclei show near-complete CG methylation loss compared to 2C nuclei — a remethylation failure not seen elsewhere.
- **Small RNA pathway activity in phloem** — Phloem companion cells express high levels of small RNA biogenesis genes (the upstream half of the RdDM pathway) but low levels of the downstream DNA methylation machinery. This suggests phloem may export small RNAs to neighboring cells rather than using them locally to maintain its own methylation.

**This repository's question:** Within this atlas, which genes show cell-type-specific gene-body CG methylation that is negatively correlated with their RNA expression? Identifying these genes connects the observed epigenomic remodelling to transcriptional regulation across the 17 cluster atlas. The analysis uses weighted least-squares (WLS) regression — methylation z-scores as the predictor, RNA expression z-scores as the outcome, weighted by sequencing coverage — run independently for each gene across the 17 clusters.

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
