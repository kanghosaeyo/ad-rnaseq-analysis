"""
Secondary differential expression comparisons within the cohort.

Runs two additional pairwise comparisons beyond the primary contrast:
one unpaired comparison against the reference group, and one paired
comparison between the two conditions measured on the same individuals.
The paired comparison blocks on patient identity to account for
between-person baseline variation, which the unpaired comparison cannot
do since its groups don't share individuals.

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
TABLES_DIR = REPO_ROOT / "results" / "tables"
FIGURES_DIR = REPO_ROOT / "results" / "figures"

EXCEL_CORRUPTED_GENE_SYMBOLS = ["1-Mar", "2-Mar"]
MIN_TOTAL_COUNT = 10
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 1.0


def load_data():
    counts = pd.read_csv(COUNTS_PATH, index_col=0).T
    meta = pd.read_csv(METADATA_PATH, index_col="sample_id")
    meta = meta.loc[counts.index]
    present = [g for g in EXCEL_CORRUPTED_GENE_SYMBOLS if g in counts.columns]
    if present:
        counts = counts.drop(columns=present)
    return counts, meta


def run_contrast(counts: pd.DataFrame, meta: pd.DataFrame, design: str,
                  contrast: list, min_count: int = MIN_TOTAL_COUNT) -> pd.DataFrame:
    keep = counts.sum(axis=0) >= min_count
    counts_f = counts.loc[:, keep]

    dds = DeseqDataSet(counts=counts_f, metadata=meta, design=design)
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=contrast)
    stat_res.summary()

    results = stat_res.results_df.copy()
    results["significant"] = (
        (results["padj"] < PADJ_THRESHOLD)
        & (results["log2FoldChange"].abs() > LFC_THRESHOLD)
    )
    return results.sort_values("padj")


def plot_volcano(results: pd.DataFrame, title: str, x_label: str, out_path: Path):
    plt.figure(figsize=(7, 6))
    not_sig = results[~results["significant"]]
    sig = results[results["significant"]]

    plt.scatter(not_sig["log2FoldChange"], -np.log10(not_sig["padj"]),
                s=8, alpha=0.3, color="gray", label="not significant")
    plt.scatter(sig["log2FoldChange"], -np.log10(sig["padj"]),
                s=10, alpha=0.6, color="crimson", label="significant")
    plt.axhline(-np.log10(PADJ_THRESHOLD), color="black", linestyle="--", linewidth=0.5)
    plt.axvline(LFC_THRESHOLD, color="black", linestyle="--", linewidth=0.5)
    plt.axvline(-LFC_THRESHOLD, color="black", linestyle="--", linewidth=0.5)

    top_genes = sig.sort_values("padj").head(8)
    for gene, row in top_genes.iterrows():
        plt.annotate(gene, (row["log2FoldChange"], -np.log10(row["padj"])),
                     fontsize=8, xytext=(3, 3), textcoords="offset points")

    plt.xlabel(x_label)
    plt.ylabel("-log10(adjusted p-value)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    print(f"Saved volcano plot -> {out_path}")
    plt.close()


def run_comparison(name: str, counts: pd.DataFrame, meta: pd.DataFrame,
                    design: str, contrast: list, plot_title: str, x_label: str):
    print(f"\n=== {name} ===")
    results = run_contrast(counts, meta, design, contrast)

    n_sig = results["significant"].sum()
    print(f"Significant genes (padj < {PADJ_THRESHOLD}, |log2FC| > {LFC_THRESHOLD}): "
          f"{n_sig} / {len(results)}")
    print(results.head(5)[["baseMean", "log2FoldChange", "padj"]])

    table_path = TABLES_DIR / f"de_{name}.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(table_path)
    print(f"Saved DE results table -> {table_path}")

    figure_path = FIGURES_DIR / f"volcano_{name}.png"
    plot_volcano(results, plot_title, x_label, figure_path)


def main():
    counts, meta = load_data()

    # Comparison 1: non-lesional vs. healthy — unpaired (different individuals),
    # fit across all samples in the cohort for maximum dispersion-estimation power.
    run_comparison(
        name="nonlesional_vs_healthy",
        counts=counts,
        meta=meta,
        design="~condition",
        contrast=["condition", "non-lesional", "healthy"],
        plot_title="Differential expression: non-lesional vs. healthy skin",
        x_label="log2 fold change (non-lesional vs. healthy)",
    )

    # Comparison 2: lesional vs. non-lesional — paired (same individuals),
    # restricted to the AD-only subset, blocking on patient identity.
    ad_mask = meta["group"] == "AD"
    counts_ad = counts.loc[ad_mask]
    meta_ad = meta.loc[ad_mask]
    run_comparison(
        name="lesional_vs_nonlesional",
        counts=counts_ad,
        meta=meta_ad,
        design="~patient_id + condition",
        contrast=["condition", "lesional", "non-lesional"],
        plot_title="Differential expression: lesional vs. non-lesional skin (paired)",
        x_label="log2 fold change (lesional vs. non-lesional)",
    )


if __name__ == "__main__":
    main()