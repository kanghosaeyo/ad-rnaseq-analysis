"""
Primary differential expression analysis between two conditions in the
cohort.

Fits a statistical model across all samples to estimate gene-level
variability, then extracts the specific pairwise comparison of interest.
Produces a ranked results table along with volcano and MA plots to
visualize effect size and significance across all tested genes.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

COUNTS_PATH = REPO_ROOT / "data" / "processed" / "counts_subset.csv"
METADATA_PATH = REPO_ROOT / "data" / "processed" / "metadata.csv"

OUT_DE_TABLE = REPO_ROOT / "results" / "tables" / "de_lesional_vs_healthy.csv"
OUT_VOLCANO = REPO_ROOT / "results" / "figures" / "volcano_lesional_vs_healthy.png"
OUT_MA = REPO_ROOT / "results" / "figures" / "ma_lesional_vs_healthy.png"

EXCEL_CORRUPTED_GENE_SYMBOLS = ["1-Mar", "2-Mar"]
MIN_TOTAL_COUNT = 10
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 1.0  # |log2FoldChange| > 1  ==  >2x change in either direction


def load_data():
    counts = pd.read_csv(COUNTS_PATH, index_col=0).T  # samples x genes
    meta = pd.read_csv(METADATA_PATH, index_col="sample_id")
    meta = meta.loc[counts.index]
    return counts, meta


def clean_counts(counts: pd.DataFrame) -> pd.DataFrame:
    present = [g for g in EXCEL_CORRUPTED_GENE_SYMBOLS if g in counts.columns]
    if present:
        counts = counts.drop(columns=present)
    keep = counts.sum(axis=0) >= MIN_TOTAL_COUNT
    return counts.loc[:, keep]


def run_de(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition")
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=["condition", "lesional", "healthy"])
    stat_res.summary()

    results = stat_res.results_df.copy()
    results["significant"] = (
        (results["padj"] < PADJ_THRESHOLD)
        & (results["log2FoldChange"].abs() > LFC_THRESHOLD)
    )
    return results.sort_values("padj")


def plot_volcano(results: pd.DataFrame):
    plt.figure(figsize=(7, 6))
    not_sig = results[~results["significant"]]
    sig = results[results["significant"]]

    plt.scatter(
        not_sig["log2FoldChange"], -np.log10(not_sig["padj"]),
        s=8, alpha=0.3, color="gray", label="not significant",
    )
    plt.scatter(
        sig["log2FoldChange"], -np.log10(sig["padj"]),
        s=10, alpha=0.6, color="crimson", label="significant",
    )
    plt.axhline(-np.log10(PADJ_THRESHOLD), color="black", linestyle="--", linewidth=0.5)
    plt.axvline(LFC_THRESHOLD, color="black", linestyle="--", linewidth=0.5)
    plt.axvline(-LFC_THRESHOLD, color="black", linestyle="--", linewidth=0.5)

    top_genes = sig.sort_values("padj").head(8)
    for gene, row in top_genes.iterrows():
        plt.annotate(
            gene, (row["log2FoldChange"], -np.log10(row["padj"])),
            fontsize=8, xytext=(3, 3), textcoords="offset points",
        )

    plt.xlabel("log2 fold change (lesional vs. healthy)")
    plt.ylabel("-log10(adjusted p-value)")
    plt.title("Differential expression: AD lesional vs. healthy skin")
    plt.legend()
    plt.tight_layout()
    OUT_VOLCANO.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_VOLCANO, dpi=150)
    print(f"Saved volcano plot -> {OUT_VOLCANO}")
    plt.close()


def plot_ma(results: pd.DataFrame):
    plt.figure(figsize=(7, 6))
    not_sig = results[~results["significant"]]
    sig = results[results["significant"]]

    plt.scatter(
        np.log10(not_sig["baseMean"] + 1), not_sig["log2FoldChange"],
        s=8, alpha=0.3, color="gray", label="not significant",
    )
    plt.scatter(
        np.log10(sig["baseMean"] + 1), sig["log2FoldChange"],
        s=10, alpha=0.6, color="crimson", label="significant",
    )
    plt.axhline(0, color="black", linewidth=0.5)
    plt.xlabel("log10(mean expression + 1)")
    plt.ylabel("log2 fold change (lesional vs. healthy)")
    plt.title("MA plot: AD lesional vs. healthy skin")
    plt.legend()
    plt.tight_layout()
    OUT_MA.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_MA, dpi=150)
    print(f"Saved MA plot -> {OUT_MA}")
    plt.close()


def main():
    counts, meta = load_data()
    counts = clean_counts(counts)
    print(f"Running DE on {counts.shape[0]} samples x {counts.shape[1]} genes")

    results = run_de(counts, meta)

    n_sig = results["significant"].sum()
    print(f"\nSignificant genes (padj < {PADJ_THRESHOLD}, "
          f"|log2FC| > {LFC_THRESHOLD}): {n_sig} / {len(results)}")
    print("\nTop 10 genes by adjusted p-value:")
    print(results.head(10)[["baseMean", "log2FoldChange", "padj"]])

    OUT_DE_TABLE.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT_DE_TABLE)
    print(f"\nSaved DE results table -> {OUT_DE_TABLE}")

    plot_volcano(results)
    plot_ma(results)


if __name__ == "__main__":
    main()