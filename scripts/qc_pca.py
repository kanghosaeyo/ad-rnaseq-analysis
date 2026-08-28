"""
Quality control and exploratory analysis on the curated cohort.

Checks sequencing depth per sample, removes any gene identifiers known to
be unreliable, filters out low-signal genes, and normalizes the count
matrix. Runs PCA to visualize overall sample similarity and flags samples
that fall unusually far from their expected group, so they can be
reviewed before downstream statistical testing.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from pydeseq2.dds import DeseqDataSet

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

COUNTS_PATH = REPO_ROOT / "data" / "processed" / "counts_subset.csv"
METADATA_PATH = REPO_ROOT / "data" / "processed" / "metadata.csv"

OUT_LIBSIZE_TABLE = REPO_ROOT / "results" / "tables" / "library_size_qc.csv"
OUT_PCA_TABLE = REPO_ROOT / "results" / "tables" / "pca_coordinates.csv"
OUT_PCA_FIGURE = REPO_ROOT / "results" / "figures" / "pca_plot.png"

# Gene symbols known to be corrupted by spreadsheet date auto-formatting
EXCEL_CORRUPTED_GENE_SYMBOLS = ["1-Mar", "2-Mar"]

MIN_TOTAL_COUNT = 10          # gene-level filter: min summed counts across samples
LOW_LIBRARY_FRACTION = 0.3    # flag samples below this fraction of the median library size
PCA_OUTLIER_SD = 2.0          # flag samples beyond this many SDs from their group centroid


def load_data():
    counts = pd.read_csv(COUNTS_PATH, index_col=0).T  # samples x genes
    meta = pd.read_csv(METADATA_PATH, index_col="sample_id")
    meta = meta.loc[counts.index]
    return counts, meta


def check_library_sizes(counts: pd.DataFrame) -> pd.DataFrame:
    lib_sizes = counts.sum(axis=1)
    median = lib_sizes.median()
    qc = pd.DataFrame({
        "sample_id": lib_sizes.index,
        "library_size": lib_sizes.values,
        "fraction_of_median": lib_sizes.values / median,
    })
    low = qc[qc["fraction_of_median"] < LOW_LIBRARY_FRACTION]
    if len(low) > 0:
        print(f"WARNING: {len(low)} sample(s) below "
              f"{LOW_LIBRARY_FRACTION:.0%} of median library size:")
        print(low["sample_id"].tolist())
    else:
        print(f"Library sizes OK: min={lib_sizes.min():,.0f}, "
              f"median={median:,.0f}, max={lib_sizes.max():,.0f}. "
              f"No samples below {LOW_LIBRARY_FRACTION:.0%} of median.")
    return qc


def drop_corrupted_genes(counts: pd.DataFrame) -> pd.DataFrame:
    present = [g for g in EXCEL_CORRUPTED_GENE_SYMBOLS if g in counts.columns]
    if present:
        print(f"Dropping {len(present)} Excel-corrupted gene symbol(s): {present}")
        counts = counts.drop(columns=present)
    return counts


def filter_low_count_genes(counts: pd.DataFrame) -> pd.DataFrame:
    keep = counts.sum(axis=0) >= MIN_TOTAL_COUNT
    print(f"Gene filter: keeping {keep.sum():,} / {len(keep):,} genes "
          f"(total count >= {MIN_TOTAL_COUNT})")
    return counts.loc[:, keep]


def normalize_and_vst(counts: pd.DataFrame, meta: pd.DataFrame) -> np.ndarray:
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition")
    dds.fit_size_factors()
    dds.vst()
    return dds.layers["vst_counts"]


def run_pca(vst_counts: np.ndarray, meta: pd.DataFrame) -> pd.DataFrame:
    pca = PCA(n_components=10)
    coords = pca.fit_transform(vst_counts)
    pca_df = pd.DataFrame(
        coords[:, :2], columns=["PC1", "PC2"], index=meta.index
    )
    pca_df["condition"] = meta["condition"].values
    pca_df["patient_id"] = meta["patient_id"].values
    var_explained = pca.explained_variance_ratio_[:2]
    print(f"PC1 explains {var_explained[0]:.1%} of variance, "
          f"PC2 explains {var_explained[1]:.1%}")
    return pca_df, var_explained


def flag_pca_outliers(pca_df: pd.DataFrame):
    for condition, group in pca_df.groupby("condition"):
        centroid = group[["PC1", "PC2"]].mean()
        dists = np.sqrt(((group[["PC1", "PC2"]] - centroid) ** 2).sum(axis=1))
        threshold = dists.mean() + PCA_OUTLIER_SD * dists.std()
        outliers = group.index[dists > threshold].tolist()
        if outliers:
            print(f"WARNING: potential PCA outlier(s) in '{condition}' group: {outliers}")


def plot_pca(pca_df: pd.DataFrame, var_explained: np.ndarray):
    plt.figure(figsize=(7, 6))
    sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="condition",
        style="condition", s=80, alpha=0.85,
    )
    plt.xlabel(f"PC1 ({var_explained[0]:.1%} variance)")
    plt.ylabel(f"PC2 ({var_explained[1]:.1%} variance)")
    plt.title("PCA of VST-normalized expression: AD lesional / non-lesional / healthy")
    plt.tight_layout()
    OUT_PCA_FIGURE.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_PCA_FIGURE, dpi=150)
    print(f"Saved PCA plot -> {OUT_PCA_FIGURE}")


def main():
    counts, meta = load_data()
    print(f"Loaded {counts.shape[0]} samples x {counts.shape[1]} genes")

    libsize_qc = check_library_sizes(counts)
    OUT_LIBSIZE_TABLE.parent.mkdir(parents=True, exist_ok=True)
    libsize_qc.to_csv(OUT_LIBSIZE_TABLE, index=False)
    print(f"Saved library size QC table -> {OUT_LIBSIZE_TABLE}")

    counts = drop_corrupted_genes(counts)
    counts = filter_low_count_genes(counts)

    vst_counts = normalize_and_vst(counts, meta)
    pca_df, var_explained = run_pca(vst_counts, meta)

    flag_pca_outliers(pca_df)

    pca_df.to_csv(OUT_PCA_TABLE)
    print(f"Saved PCA coordinates -> {OUT_PCA_TABLE}")

    plot_pca(pca_df, var_explained)


if __name__ == "__main__":
    main()