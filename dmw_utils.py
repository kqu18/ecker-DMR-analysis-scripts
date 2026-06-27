"""Shared helpers for the methylation × RNA pipeline.

Pipeline stages this module covers (see CLAUDE.md for the full project context):

    1. Load RNA + cluster metadata        -> load_rna_cluster_matrix
    2. Load gene or DMW MCDS              -> load_gene_mcds / load_dmw_mcds
    3. Per-cluster mc/cov aggregation     -> aggregate_per_cluster
    4. Logit-adjusted methylation residual -> compute_logit_residuals
    5. DMW -> gene-body overlap            -> map_dmw_to_genes
    6. WLS correlation across clusters     -> wls_correlate
    7. Top-N anticorrelation heatmap       -> plot_paired_heatmap
    8. Cluster-specific hypomethylation    -> identify_relative_lows
    9. Multi-genotype residual matrices    -> derive_genotype_from_cellid,
                                              compute_matrices_by_genotype,
                                              plot_genotype_panels

The functions are factored out of the canonical notebooks:
    notebooks/2026/260305_corr_heatmap_refactored.ipynb   (gene-level WLS)
    notebooks/2026/260325_identify_hypo.ipynb             (gene hypo)
    notebooks/2026/260428_dmw_corr_heatmap_claude.ipynb   (DMW-level WLS)
    notebooks/2026/260501_dmw_identify_hypo.ipynb         (DMW hypo)
    notebooks/2026/260514_dmw_methylation_matrices_by_genotype.ipynb
"""

from __future__ import annotations

import os
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
import statsmodels.api as sm
from matplotlib.colors import LinearSegmentedColormap
from scipy.special import logit, expit
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests


# ============================================================================
# Constants
# ============================================================================

CLUSTER_MAPPING: dict = {
    3: '3 Early M', 0: '0 Mid M', 1: '1 Late M', 7: '7 Expanding M', 10: '10 Senescent M',
    5: '5 Early E', 4: '4 Abaxial E', 6: '6 Adaxial E', 13: '13 Expanding E', 8: '8 Senescent E',
    14: '14 Guard cell', 2: '2 Vasculature', 11: '11 PPP', 9: '9 Phloem',
    16: '16 Myrosinase', 12: '12 S phase', 15: '15 G2M phase',
}
CLUSTER_MAPPING.update({str(k): v for k, v in CLUSTER_MAPPING.items()})

CLUSTER_ORDER: list = [
    '3 Early M', '0 Mid M', '7 Expanding M', '1 Late M', '10 Senescent M',
    '5 Early E', '13 Expanding E', '4 Abaxial E', '6 Adaxial E', '8 Senescent E',
    '14 Guard cell', '2 Vasculature', '11 PPP', '9 Phloem',
    '16 Myrosinase', '12 S phase', '15 G2M phase',
]

# DMW MCDS encodes genotype in cell-ID prefix `240614_mct_<G>_<plate>_<well>`
GENOTYPE_MAP: dict = {"1": "col", "2": "rdd", "3": "met"}
GENOTYPES: list = ["col", "rdd", "met"]

# Canonical paths (override per-call if needed)
PATHS = {
    "gene_mcds":  "/ceph/MethDev/pbio/hazel/gene_mcds/gene",
    "dmw_mcds":  ["/ceph/MethDev/pbio/andy/JW/section_4/All_mcds_dmw3_0.mcds",
                  "/ceph/MethDev/pbio/andy/JW/section_4/All_mcds_dmw3_1.mcds"],
    "rna":       "/ceph/MethDev/pbio/hazel/count_Arab13_with_TEs.csv",
    "meta":      "/ceph/MethDev/pbio/hazel/merged_cluster_assignments.csv",
    "genes_bed": "/ceph/MethDev/pbio/andy/JW/section_4/genes.fixed.bed",
    "chunks_gff":"/ceph/MethDev/pbio/andy/JW/section_4/chunks_CG_minfilt.gff",
    "fig_outdir":"/ceph/MethDev/pbio/kay/figures/heatmaps",
    "data_outdir":"/ceph/MethDev/pbio/kay/data",
}


# ============================================================================
# Stage 1-2 — Data loading
# ============================================================================

def load_rna_cluster_matrix(rna_path: str = PATHS["rna"],
                            meta_path: str = PATHS["meta"]):
    """RNA counts CSV (genes × cells, transposed on load) -> log1p-CPM mean per cluster.

    Returns
    -------
    cluster_rna_matrix : DataFrame (genes × clusters)
    meta              : DataFrame indexed by cell_id with `merged_cluster` column
    """
    rna  = pd.read_csv(rna_path).T
    meta = pd.read_csv(meta_path).set_index("cell_id")

    common = rna.index.intersection(meta.index)
    rna    = rna.loc[common]
    rna    = rna.div(rna.sum(axis=1), axis=0) * 1e6   # CPM
    rna    = np.log1p(rna)
    rna["cluster"] = meta.loc[common, "merged_cluster"]
    return rna.groupby("cluster").mean().T, meta


def load_gene_mcds(mcds_path: str = PATHS["gene_mcds"]):
    """Open the gene-level MCDS (var_dim='gene')."""
    from ALLCools.mcds import MCDS
    return MCDS.open(mcds_paths=mcds_path, var_dim="gene")


def load_dmw_mcds(mcds_paths: Iterable[str] = PATHS["dmw_mcds"],
                  curated_chunks_gff: Optional[str] = None):
    """Open both DMW MCDS files; optionally restrict to a curated chunk set.

    Parameters
    ----------
    mcds_paths : list of two .mcds files; ALLCools concatenates on cell.
    curated_chunks_gff : path to chunks_CG_minfilt.gff (0-based, despite .gff name).
        If given, MCDS is subset to those DMWs.
    """
    from ALLCools.mcds import MCDS
    mcds = MCDS.open(mcds_paths=list(mcds_paths), var_dim="dmw")

    if curated_chunks_gff is not None:
        chunks = pd.read_csv(
            curated_chunks_gff, sep="\t", header=None, comment="#",
            usecols=[0, 3, 4], names=["chrom", "start", "end"],
        )
        chunks["chrom"] = chunks["chrom"].str.replace("^chr", "", regex=True)
        chunks["start"] = chunks["start"].astype(int) + 1    # 0-based -> 1-based MCDS
        chunks["end"]   = chunks["end"].astype(int)

        mcds_bed = pd.DataFrame({
            "chrom":  mcds["dmw_chrom"].values.astype(str),
            "start":  mcds["dmw_start"].values.astype(int),
            "end":    mcds["dmw_end"].values.astype(int),
            "dmw_id": mcds.get_index("dmw").values,
        })
        keep = mcds_bed.merge(chunks, on=["chrom", "start", "end"], how="inner")
        mcds = mcds.sel(dmw=keep["dmw_id"].tolist())

    return mcds


def derive_genotype_from_cellid(cell_ids: Iterable[str]) -> pd.Series:
    """Map `240614_mct_<G>_<plate>_<well>` -> {col, rdd, met}. Index = cell_ids."""
    cell_ids = list(cell_ids)
    return (pd.Series(cell_ids, index=cell_ids)
              .str.split("_").str[2]
              .map(GENOTYPE_MAP))


# ============================================================================
# Stage 3-4 — Per-cluster aggregation + logit residuals
# ============================================================================

def aggregate_per_cluster(mcds, meta: pd.DataFrame,
                          var_dim: str = "gene",
                          mc_type: str = "CGN",
                          cell_subset: Optional[Iterable[str]] = None):
    """Sum mc and cov per (var_dim, merged_cluster).

    Returns (cluster_c, cluster_m), both shaped (var_dim × clusters).
    """
    da_name  = f"{var_dim}_da"
    mc_cell  = mcds.sel(mc_type=mc_type, count_type="mc")[da_name].to_pandas()
    cov_cell = mcds.sel(mc_type=mc_type, count_type="cov")[da_name].to_pandas()

    if cell_subset is not None:
        cell_subset = list(cell_subset)
        mc_cell  = mc_cell.loc[mc_cell.index.intersection(cell_subset)]
        cov_cell = cov_cell.loc[cov_cell.index.intersection(cell_subset)]

    keep = mc_cell.index.intersection(meta.index)
    mc_cell, cov_cell = mc_cell.loc[keep], cov_cell.loc[keep]
    mc_cell["cluster"]  = meta.loc[keep, "merged_cluster"]
    cov_cell["cluster"] = meta.loc[keep, "merged_cluster"]

    cluster_c = mc_cell.groupby("cluster").sum().T
    cluster_m = cov_cell.groupby("cluster").sum().T
    return cluster_c, cluster_m


def compute_logit_residuals(cluster_c: pd.DataFrame,
                            cluster_m: pd.DataFrame,
                            min_cov: int = 50,
                            high_meth_th: float = 0.10):
    """Logit-adjusted residuals: per-cluster offsets + per-item baseline.

    Filters in this order:
        1. Per-item summed coverage >= min_cov.
        2. (optional, for filtered_meth_df only) raw fraction >= high_meth_th in any cluster.

    Returns
    -------
    obs_p             : observed mc/cov per (item, cluster), coverage-filtered
    adj_meth_matrix   : obs_p - p0(cluster-adjusted)
    filtered_meth_df  : adj_meth_matrix restricted to high-methylation items
    deltas            : per-cluster logit offsets (Series, indexed by cluster)
    """
    global_c = cluster_c.sum(axis=0)
    global_m = cluster_m.sum(axis=0)
    global_p = global_c / global_m
    d        = logit(np.clip(global_p, 1e-6, 1 - 1e-6))
    deltas   = d - np.average(d, weights=global_m)

    item_m = cluster_m.sum(axis=1)
    valid  = item_m >= min_cov
    cluster_c = cluster_c[valid]
    cluster_m = cluster_m[valid]

    pbar  = np.clip(cluster_c.sum(axis=1) / cluster_m.sum(axis=1), 1e-6, 1 - 1e-6)
    p0    = pd.DataFrame(
        expit(logit(pbar).values[:, None] + deltas.values[None, :]),
        index=cluster_c.index, columns=cluster_c.columns,
    )
    obs_p = (cluster_c / cluster_m.replace(0, np.nan)).fillna(0)
    adj   = obs_p - p0

    highly_meth = obs_p[(obs_p >= high_meth_th).any(axis=1)].index
    filtered    = adj.loc[highly_meth]
    return obs_p, adj, filtered, deltas


# ============================================================================
# Stage 5 — DMW -> gene-body overlap
# ============================================================================

def map_dmw_to_genes(mcds, genes_bed_path: str = PATHS["genes_bed"]) -> pd.DataFrame:
    """Return (dmw_id, gene) pairs for every overlap between a DMW and a gene body.

    Harmonizes the chr-prefix gap: DMW MCDS uses bare ints, gene BED uses chr1-style.
    """
    import bioframe as bf

    dmw_bed = pd.DataFrame({
        "chrom":  mcds["dmw_chrom"].values,
        "start":  mcds["dmw_start"].values.astype(int),
        "end":    mcds["dmw_end"].values.astype(int),
        "dmw_id": mcds.get_index("dmw").values,
    })
    dmw_bed["chrom"] = "chr" + dmw_bed["chrom"].astype(str)

    genes = pd.read_csv(
        genes_bed_path, sep="\t", header=None,
        names=["chrom", "start", "end", "gene", "score", "strand"],
        usecols=[0, 1, 2, 3],
    )

    overlaps = bf.overlap(
        dmw_bed, genes,
        cols1=("chrom", "start", "end"),
        cols2=("chrom", "start", "end"),
        how="inner", suffixes=("_dmw", "_gene"),
    )
    dmw_col  = next(c for c in overlaps.columns if c.startswith("dmw_id"))
    gene_col = next(c for c in overlaps.columns if c == "gene" or c.startswith("gene_"))
    return (overlaps[[dmw_col, gene_col]]
            .rename(columns={dmw_col: "dmw_id", gene_col: "gene"})
            .drop_duplicates()
            .reset_index(drop=True))


# ============================================================================
# Stage 6 — WLS correlation
# ============================================================================

def wls_correlate(meth_df: pd.DataFrame,
                  rna_df: pd.DataFrame,
                  weights_df: pd.DataFrame,
                  pairs: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Per-row WLS of zscored RNA on zscored methylation, weighted by coverage.

    Two modes:
      - **Gene-level** (pairs is None): correlate each row of meth_df with the same-named
        row of rna_df. Returned columns: gene, rho, p_value, fdr.
      - **DMW-level**  (pairs given): expects pairs[['dmw_id','gene']]; correlates
        meth_df.loc[dmw_id] against rna_df.loc[gene]. Returned columns: dmw_id, gene, rho, p_value, fdr.
    """
    rows = []
    if pairs is None:
        common = meth_df.index.intersection(rna_df.index).intersection(weights_df.index)
        for item in common:
            r = rna_df.loc[item].values
            m = meth_df.loc[item].values
            w = weights_df.loc[item].values
            if np.std(r) == 0 or np.std(m) == 0:
                continue
            try:
                fit = sm.WLS(zscore(r), sm.add_constant(zscore(m)), weights=w).fit()
                rows.append((item, fit.params[1], fit.pvalues[1]))
            except Exception:
                continue
        out = pd.DataFrame(rows, columns=["gene", "rho", "p_value"])
    else:
        for dmw_id, gene in pairs.itertuples(index=False):
            if dmw_id not in meth_df.index or gene not in rna_df.index:
                continue
            r = rna_df.loc[gene].values
            m = meth_df.loc[dmw_id].values
            w = weights_df.loc[dmw_id].values
            if np.std(r) == 0 or np.std(m) == 0:
                continue
            try:
                fit = sm.WLS(zscore(r), sm.add_constant(zscore(m)), weights=w).fit()
                rows.append((dmw_id, gene, fit.params[1], fit.pvalues[1]))
            except Exception:
                continue
        out = pd.DataFrame(rows, columns=["dmw_id", "gene", "rho", "p_value"])

    if len(out):
        out["fdr"] = multipletests(out["p_value"], method="fdr_bh")[1]
    return out


# ============================================================================
# Stage 7 — Heatmap plotting
# ============================================================================

def _prepare_zscore(df: pd.DataFrame, clip: float = 2.5) -> pd.DataFrame:
    """Row z-score, clipped to +/- `clip`."""
    return df.apply(zscore, axis=1).clip(lower=-clip, upper=clip)


def _ordered_columns(df: pd.DataFrame) -> list:
    df = df.rename(columns=CLUSTER_MAPPING)
    return [c for c in CLUSTER_ORDER if c in df.columns]


def plot_paired_heatmap(meth_df: pd.DataFrame,
                        rna_df: pd.DataFrame,
                        title_meth: str,
                        title_rna: str = "RNA Expression (log1p CPM)",
                        outpath: Optional[str] = None,
                        rna_cmap_colors=("orange", "white", "teal"),
                        ylabel: Optional[str] = None):
    """Side-by-side meth + RNA heatmap; rows hierarchically clustered on methylation.

    Both inputs must be (rows × clusters), with matching row indices (e.g. labels
    like "DMW_123 | AT5G53210" for DMW–gene pairs, or gene names for the gene-level case).

    `ylabel` overrides the default "rows (n=...)" — use it to clarify when rows
    differ across panels (e.g. "Top-300 DMW–gene pairs (per-subset WLS, n=300)").
    """
    meth_z = _prepare_zscore(meth_df).rename(columns=CLUSTER_MAPPING)
    rna_z  = _prepare_zscore(rna_df ).rename(columns=CLUSTER_MAPPING)
    cols   = _ordered_columns(meth_df)
    meth_z, rna_z = meth_z[cols], rna_z[cols]

    link  = hc.linkage(sp.distance.pdist(meth_z.values), method="ward")
    order = hc.leaves_list(link)
    meth_z, rna_z = meth_z.iloc[order], rna_z.iloc[order]

    rna_cmap = LinearSegmentedColormap.from_list("rna_cmap", list(rna_cmap_colors))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 12),
                                   gridspec_kw={"width_ratios": [1, 1]})

    sns.heatmap(meth_z, cmap="coolwarm", center=0, ax=ax1, yticklabels=False,
                cbar_kws={"label": "Methylation Residual (Z-score)"})
    sns.heatmap(rna_z, cmap=rna_cmap, ax=ax2, yticklabels=False,
                cbar_kws={"label": "RNA Expression (Z-score)"})

    ax1.set_title(title_meth, fontsize=15, pad=15)
    ax2.set_title(title_rna,  fontsize=15, pad=15)
    ax1.set_ylabel(ylabel if ylabel is not None else f"rows (n={len(meth_z)})",
                   fontsize=12)
    for ax in (ax1, ax2):
        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right",
                           rotation_mode="anchor")
        for i, col in enumerate(cols):
            if "Senescent M" in col or "Senescent E" in col:
                ax.axvline(i + 1, color="black", lw=1.2, linestyle=":")

    plt.tight_layout()
    if outpath is not None:
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
    return fig


# ============================================================================
# Stage 8 — Cluster-specific hypomethylation
# ============================================================================

def identify_relative_lows(item_list: Iterable,
                           obs_p: pd.DataFrame,
                           adj_meth_matrix: pd.DataFrame,
                           cluster_m: pd.DataFrame,
                           min_coverage: int = 10,
                           residual_threshold: float = -0.05,
                           zscore_threshold: float = -1.0,
                           min_raw_max: float = 0.05) -> pd.DataFrame:
    """Flag items (genes or DMWs) whose methylation is cluster-specifically low.

    Item qualifies in a cluster when:
      - Raw fraction exceeds `min_raw_max` in at least one cluster (globally methylated)
      - >= 3 clusters have coverage >= `min_coverage`
      - Residual <= `residual_threshold` (absolute drop)
      - Within-item z-score <= `zscore_threshold` (relative drop)
    """
    rows = []
    valid = [i for i in item_list if i in adj_meth_matrix.index]

    for item in valid:
        obs = obs_p.loc[item]
        adj = adj_meth_matrix.loc[item]
        cov = cluster_m.loc[item]

        if obs.max() < min_raw_max:
            continue
        mask = cov >= min_coverage
        if mask.sum() < 3:
            continue

        adj_v = adj[mask]
        zs    = pd.Series(zscore(adj_v), index=adj_v.index)
        hits  = adj_v[(adj_v <= residual_threshold) & (zs <= zscore_threshold)].index

        for cluster in hits:
            rows.append({
                "item": item, "cluster": cluster,
                "raw_fraction": obs[cluster], "residual": adj_v[cluster],
                "z_score": zs[cluster],       "coverage": cov[cluster],
            })

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(by=["item", "residual"])


# ============================================================================
# Stage 9 — Multi-genotype matrices and panels
# ============================================================================

def compute_matrices_by_genotype(mcds, meta: pd.DataFrame,
                                 mc_type: str = "CGN",
                                 min_cov: int = 50,
                                 high_meth_th: float = 0.10,
                                 canonical_rows_from: str = "col"):
    """Per-genotype raw + residual methylation matrices on a shared row set.

    Per-cluster deltas are RE-FIT per genotype so each residual matrix shows
    within-genotype cluster structure (not vs col-0). The canonical row set is
    derived from `canonical_rows_from` (default col-0): items >= high_meth_th
    in any cluster.

    Returns
    -------
    raw_by_geno : dict {genotype: DataFrame raw obs_p}
    res_by_geno : dict {genotype: DataFrame residual}
    canonical_rows : Index of items kept (anchored on `canonical_rows_from`)
    """
    geno_per_cell = derive_genotype_from_cellid(mcds["cell"].values)
    var_dim = "dmw" if "dmw" in mcds.dims else "gene"

    raw_by_geno, res_by_geno = {}, {}
    for g in GENOTYPES:
        cells_g = geno_per_cell.index[geno_per_cell == g]
        c, m = aggregate_per_cluster(mcds, meta, var_dim=var_dim,
                                     mc_type=mc_type, cell_subset=cells_g)
        obs_p, adj, _, _ = compute_logit_residuals(c, m,
                                                   min_cov=min_cov,
                                                   high_meth_th=0)  # no high-meth filter yet
        raw_by_geno[g] = obs_p
        res_by_geno[g] = adj

    # Canonical row set from the anchor genotype (typically col-0)
    canonical = raw_by_geno[canonical_rows_from]
    canonical_rows = canonical[(canonical >= high_meth_th).any(axis=1)].index

    for g in GENOTYPES:
        raw_by_geno[g] = raw_by_geno[g].reindex(canonical_rows)
        res_by_geno[g] = res_by_geno[g].reindex(canonical_rows)

    return raw_by_geno, res_by_geno, canonical_rows


def plot_genotype_panels(residuals_by_geno: dict,
                         row_order: Optional[pd.Index] = None,
                         vmax: float = 2.5,
                         outpath: Optional[str] = None,
                         title: str = "",
                         zscore_rows: bool = True,
                         cbar_label: Optional[str] = None,
                         show_bulk_shift: bool = False):
    """Side-by-side col / rdd / met residual heatmaps.

    Parameters
    ----------
    zscore_rows : True for per-genotype residuals (default), False for col-anchored
        residuals where the bulk shift between genotypes should be preserved.
    show_bulk_shift : append mean residual per panel to its title (use with the
        col-anchored mode so rdd's red shift and met's blue shift are quantified).

    If `row_order` is None, uses col's Ward dendrogram on z-scored residuals.
    Rows containing any NaN across the three panels are dropped for alignment.
    """
    cols = _ordered_columns(residuals_by_geno["col"])
    res  = {g: residuals_by_geno[g].rename(columns=CLUSTER_MAPPING)[cols]
            for g in GENOTYPES}

    nan_free = set(res["col"].dropna().index)
    for g in GENOTYPES:
        nan_free &= set(res[g].dropna().index)

    if row_order is None:
        col_z   = _prepare_zscore(res["col"].loc[list(nan_free)])
        link    = hc.linkage(col_z.values, method="ward", metric="euclidean")
        row_order = col_z.index[hc.leaves_list(link)]
    row_order = pd.Index([r for r in row_order if r in nan_free])

    if cbar_label is None:
        cbar_label = ("Methylation residual (row Z-score)" if zscore_rows
                      else "obs_p − p0(col)")

    fig, axes = plt.subplots(1, 3, figsize=(28, 14), sharey=True,
                             gridspec_kw={"wspace": 0.05})
    for ax, g in zip(axes, GENOTYPES):
        mat = res[g].loc[row_order]
        if zscore_rows:
            mat = _prepare_zscore(mat)
        sns.heatmap(mat, ax=ax, cmap="coolwarm", center=0, vmin=-vmax, vmax=vmax,
                    yticklabels=False, cbar=(g == GENOTYPES[-1]),
                    cbar_kws={"label": cbar_label})
        panel_title = g
        if show_bulk_shift:
            panel_title += f"   (mean shift vs col-0: {float(np.nanmean(mat.values)):+.3f})"
        ax.set_title(panel_title, fontsize=14)
        ax.set_xlabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right",
                           rotation_mode="anchor")
        for i, col in enumerate(mat.columns):
            if "Senescent M" in col or "Senescent E" in col:
                ax.axvline(i + 1, color="black", lw=0.8, linestyle=":")

    axes[0].set_ylabel(f"rows (n={len(row_order)})", fontsize=12)
    if title:
        fig.suptitle(title, fontsize=15, y=1.02)
    if outpath is not None:
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        fig.savefig(outpath, dpi=300, bbox_inches="tight")
    return fig, row_order


def compute_col_anchored_residuals(mcds, meta: pd.DataFrame,
                                   canonical_rows: pd.Index,
                                   mc_type: str = "CGN",
                                   min_cov: int = 50) -> dict:
    """For each genotype, return obs_p[geno] - p0(col).

    Per-cluster deltas and per-DMW pbar are fit ONCE on col-0, then every genotype
    is scored against that same null. So:
      - col panel sits near zero (matches its own null)
      - rdd panel shifts red (RdDM mutants are hypermethylated vs col)
      - met panel shifts blue (MET1 mutants are hypomethylated vs col)

    Output rows are aligned to `canonical_rows` (typically the col-0 high-meth set
    returned by `compute_matrices_by_genotype`). Rows with insufficient coverage
    in a given genotype get NaN.
    """
    geno_per_cell = derive_genotype_from_cellid(mcds["cell"].values)
    var_dim = "dmw" if "dmw" in mcds.dims else "gene"

    col_cells = geno_per_cell.index[geno_per_cell == "col"]
    col_c, col_m = aggregate_per_cluster(mcds, meta, var_dim=var_dim,
                                         mc_type=mc_type, cell_subset=col_cells)

    g_p     = col_c.sum(axis=0) / col_m.sum(axis=0)
    d       = logit(np.clip(g_p, 1e-6, 1 - 1e-6))
    deltas  = d - np.average(d, weights=col_m.sum(axis=0))
    pbar    = np.clip(col_c.sum(axis=1) / col_m.sum(axis=1), 1e-6, 1 - 1e-6)
    col_p0  = pd.DataFrame(
        expit(logit(pbar).values[:, None] + deltas.values[None, :]),
        index=col_c.index, columns=col_c.columns,
    )

    anchored = {}
    for g in GENOTYPES:
        cells_g = geno_per_cell.index[geno_per_cell == g]
        c, m = aggregate_per_cluster(mcds, meta, var_dim=var_dim,
                                     mc_type=mc_type, cell_subset=cells_g)
        obs_p = c / m.replace(0, np.nan)
        obs_p = obs_p.where(m.sum(axis=1) >= min_cov, np.nan)
        anchored[g] = (obs_p - col_p0).reindex(canonical_rows)
    return anchored


# ============================================================================
# Cluster subsets — for the 9-condition supplementary atlas
# ============================================================================

def make_cluster_subsets(available_clusters: Iterable[str]) -> dict:
    """Return the 9 biologically-motivated cluster subsets, intersected with what's
    actually present in the data. Used to drive the supplementary atlas (Fig S1).

    Keys are stable filename-safe slugs used by 260318_9conditions.ipynb.
    """
    avail = [c for c in CLUSTER_ORDER if c in set(available_clusters)]
    meso  = [c for c in ['3 Early M', '0 Mid M', '7 Expanding M',
                         '1 Late M', '10 Senescent M'] if c in avail]
    epi   = [c for c in ['5 Early E', '13 Expanding E', '4 Abaxial E',
                         '6 Adaxial E', '8 Senescent E'] if c in avail]
    other_all  = ['14 Guard cell', '2 Vasculature', '11 PPP', '9 Phloem',
                  '16 Myrosinase', '12 S phase', '15 G2M phase']
    terminal_other = ['14 Guard cell', '2 Vasculature', '11 PPP', '9 Phloem',
                      '16 Myrosinase']
    terminal_ex_guard = ['2 Vasculature', '11 PPP', '9 Phloem', '16 Myrosinase']

    return {
        "1_All_Clusters":          avail,
        "2_Meta_States":           avail,    # plotted with meso/epi collapsed downstream
        "3_Mesophyll_Only":        meso,
        "4_Epidermis_Only":        epi,
        "5_Non_Meso_Epi":          [c for c in avail if c in other_all],
        "6_Terminal_Non_Meso_Epi": [c for c in avail if c in terminal_other],
        "7_Terminal_Ex_Guard":     [c for c in avail if c in terminal_ex_guard],
        "8_All_Ex_Guard_G2M":      [c for c in avail
                                    if c not in ('14 Guard cell', '15 G2M phase')],
        "9_Meso_and_Epi":          meso + epi,
    }
