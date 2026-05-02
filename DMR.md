# **Python Pipeline for Differential Methylation Analysis**

This document outlines a bioinformatics pipeline in Python for identifying Differentially Methylated Regions (DMRs) from clustered bisulfite sequencing data. The core of the analysis is a **Pearson's chi-square ($X^2$) test** performed on genomic windows to find regions with significant methylation variance between cell clusters.

-----

## **1. Setup: Import Libraries**

First, we import the necessary Python libraries for data manipulation, numerical operations, statistical calculations, and plotting.

```python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.special import logit, expit
from scipy.stats import chi2
from joblib import Parallel, delayed
```

-----

## **2. Data Loading and Preprocessing**

This section covers loading the annotated methylation data and performing an initial cleaning step.

### **Load Annotated Methylation Data**

The data is loaded from a tab-separated file into a pandas **DataFrame**. Each row represents a specific genomic location within a cell cluster.

  - **`c`**: Number of methylated reads (Cytosines).
  - **`t`**: Number of unmethylated reads (Thymines).
  - **`score`**: Methylation ratio ($c / (c+t)$).
  - **`cluster`**: The cell cluster ID.

<!-- end list -->

```python
# Load the dataset from a TSV file
df = pd.read_csv('../data/annotated_filtered_col.CG_2.fast.tsv', sep='\t')
```

### **Quality Control (QC)**

We filter the data by removing any rows that have a `NaN` value in the `score_masked` column. This step ensures that we only work with high-quality data points that have passed previous filtering stages.

```python
# Drop rows where 'score_masked' is NaN
df = df.dropna(subset=['score_masked'])
```

-----

## **3. Statistical Modeling and Functions**

Here we define the core functions for our statistical analysis. The main goal is to calculate a chi-square statistic for each genomic window to test for differential methylation.

### **Calculate Global Cluster Offsets (`fit_offsets`)**

This function calculates a global, baseline methylation level for each cell **cluster**. By aggregating counts across the entire genome, we can determine each cluster's typical methylation propensity. These "offsets" are calculated in logit space and are used to establish the null hypothesis for our statistical test.

```python
def fit_offsets(df_all):
    """Calculates global methylation offsets for each cluster."""
    # Aggregate methylated (c) and unmethylated (t) counts for each cluster
    agg = df_all.groupby('cluster')[['c','t']].sum()
    M = agg.sum(axis=1)
    p = agg['c'] / M

    # Convert to logit space and center the offsets
    d = logit(np.clip(p, 1e-6, 1-1e-6))
    return d - np.average(d, weights=M)
```

### **Calculate Window-wise Statistics (`window_stats`)**

This function iterates through each unique genomic window (defined by `chr`, `start`, `end`) and performs the chi-square test.

The key steps for each window are:

1.  **Filter Clusters**: Remove clusters with low coverage (fewer than 5 reads).
2.  **Establish Null Hypothesis ($p_0$)**: Calculate the expected methylation level for each cluster in the window. This is done by taking the window's overall methylation average and adjusting it by each cluster's global offset (`deltas`).
3.  **Calculate $X^2$ Statistic**: A Pearson's chi-square statistic is calculated to compare the *observed* methylated counts (`c`) to the *expected* counts under the null hypothesis (`E = m * p0`).
4.  **Calculate P-value**: The $X^2$ value is used to calculate a p-value, which is adjusted for overdispersion using a calculated `phi` factor.

<!-- end list -->

```python
def window_stats(df_all, deltas, tau=20.0, eps=1e-6):
    """Computes chi-square statistics for each genomic window."""
    out = []
    for (chr_, start, end), g in df_all.groupby(['chr','start','end']):
        c = g['c'].to_numpy(); t = g['t'].to_numpy(); m = c + t
        k = g['cluster'].to_numpy()

        # Keep clusters with sufficient coverage (>= 5 reads)
        keep = m >= 5
        if keep.sum() < 2:
            continue
        c, m, k = c[keep], m[keep], k[keep]

        # Calculate the expected methylation rate (p0) under the null hypothesis
        pbar = c.sum() / m.sum()
        p0 = expit(logit(np.clip(pbar, 1e-6, 1-1e-6)) + deltas.loc[k].to_numpy())

        # Pearson X^2 with binomial variance, epsilon-stabilized
        E = m * p0
        Var = m * p0 * (1 - p0) + eps
        X2 = ((c - E)**2 / Var).sum()
        dfree = len(c) - 1

        # Calculate max deviation and identify high/low clusters
        p_hat = c / m
        dev = np.abs(p_hat - p0)
        dmax = float(dev.max() - dev.min())
        hi = int(k[p0.argmax()])
        lo = int(k[p0.argmin()])

        out.append((chr_, start, end, X2, m.sum(), dfree, dmax, hi, lo))

    res = pd.DataFrame(out, columns=['chr','start','end','X2', '∑m', 'df','delta_max','hi_cluster','lo_cluster'])

    # Correct for overdispersion using a global phi factor
    phi = np.median(res['X2'] / np.maximum(res['df'], 1)) if len(res) else 1.0
    res['pval'] = 1 - chi2.cdf(res['X2'] / max(phi, 1e-6), res['df'])
    res['phi'] = phi
    return res
```

### **Multiple Hypothesis Correction (`bh_fdr`)**

Since we are performing thousands of tests (one for each genomic window), we must correct for multiple hypotheses. This function implements the **Benjamini-Hochberg procedure** to control the False Discovery Rate (FDR) and calculate q-values.

```python
def bh_fdr(p):
    """Benjamini-Hochberg False Discovery Rate correction."""
    if len(p) == 0:
        return p
    r = np.argsort(p)
    ranks = np.empty_like(r); ranks[r] = np.arange(1, len(p) + 1)
    q = p * len(p) / np.maximum(ranks, 1)
    q_sorted = np.minimum.accumulate(np.sort(q)[::-1])[::-1]
    out = np.empty_like(q_sorted); out[r] = q_sorted
    return np.clip(out, 0, 1)
```

-----

## **4. Running the Analysis Pipeline**

Now, we execute the functions in order to produce the final statistics table.

```python
# 1. Calculate global cluster offsets from the entire dataset
deltas = fit_offsets(df)

# 2. Compute window-by-window statistics
stats = window_stats(df, deltas, tau=20.0)

# 3. Apply FDR correction to get q-values
stats['qval'] = bh_fdr(stats['pval'].to_numpy())
```

The resulting `stats` DataFrame contains the test results for each genomic window, including the $X^2$ statistic, p-value, and q-value.

-----

## **5. Results and Visualization**

A simple way to explore the results is to visualize the distribution of the effect size, `delta_max`. This value represents the maximum difference in methylation deviation within a window, giving an idea of how strong the differential methylation is.

```python
# Plot a histogram of the maximum methylation deviation (delta_max)
plt.hist(stats['delta_max'], bins=100)
plt.xlabel("Delta Max")
plt.ylabel("Frequency")
plt.title("Distribution of Maximum Methylation Deviation per Window")
plt.show()
```

-----

## **6. Saving Results**

Finally, the resulting statistics table is saved to a CSV file for further analysis, filtering, and visualization.

```python
# Save the final statistics table to a CSV file
stats.to_csv("../data/x2_CG_stat_all_std.csv", index=False)
```

-----

## **Appendix: Parallel Processing**

For very large datasets, the window statistics calculation can be parallelized to improve performance. The function below uses `joblib` to distribute the computation across multiple CPU cores.

*Note: The notebook showed this parallel version was interrupted, suggesting the single-threaded `window_stats` function was used for the final run. However, this parallel implementation is a useful alternative for performance optimization.*

```python
def compute_window_stats(chr_, start, end, g, deltas, a0, b0, eps):
    """Helper function for a single window (for parallel execution)."""
    c = g['c'].to_numpy(); t = g['t'].to_numpy(); m = c + t
    k = g['cluster'].to_numpy()
    keep = m >= 5
    if keep.sum() < 2:
        return None
    c, m, k = c[keep], m[keep], k[keep]
    p_hat = c / m
    pbar = c.sum() / m.sum()
    p0 = expit(logit(np.clip(pbar, 1e-6, 1-1e-6)) + deltas.loc[k].to_numpy())
    dev = np.abs(p_hat - p0)
    E = m * p0
    Var = m * p0 * (1 - p0) + eps
    X2 = ((c - E)**2 / Var).sum()
    dfree = len(c) - 1
    dmax = float(dev.max() - dev.min())
    hi = int(k[p0.argmax()])
    lo = int(k[p0.argmin()])
    return (chr_, start, end, X2, m.sum(), dfree, dmax, hi, lo)

def window_stats_parallel(df_all, deltas, tau=20.0, eps=1e-6, n_jobs=8):
    """Computes chi-square statistics in parallel."""
    agg = df_all.groupby('cluster')[['c','t']].sum()
    mu = agg['c'].sum() / agg.sum(axis=1).sum()
    a0, b0 = mu * tau, (1 - mu) * tau

    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(compute_window_stats)(chr_, start, end, g, deltas, a0, b0, eps)
        for (chr_, start, end), g in df_all.groupby(['chr','start','end'])
    )

    out = [r for r in results if r is not None]
    res = pd.DataFrame(out, columns=['chr','start','end','X2','∑m','df','delta_max','hi_cluster','lo_cluster'])
    phi = np.median(res['X2'] / np.maximum(res['df'], 1)) if len(res) else 1.0
    res['pval'] = 1 - chi2.cdf(res['X2'] / max(phi, 1e-6), res['df'])
    res['phi'] = phi
    return res
```
