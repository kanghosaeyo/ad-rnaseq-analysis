"""
Targeted visualization of curated genes across the three conditions.

Builds a heatmap combining the top differentially expressed genes with a
small set of named barrier and immune genes of specific interest, plus
boxplots for a handful of standout genes. Uses VST-normalized expression
(same transform as the QC/PCA step) rather than raw counts, since raw
counts aren't directly comparable across samples with different
sequencing depth.

"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

COUNTS_PATH = REPO_ROOT / "data" / "processed" / "counts_subset.csv"
METADATA_PATH = REPO_ROOT / "data" / "processed" / "metadata.csv"
DE_PATH = REPO_ROOT / "results" / "tables" / "de_lesional_vs_healthy.csv"

OUT_HEATMAP = REPO_ROOT / "results" / "figures" / "heatmap_curated_panel.png"
OUT_BOXPLOTS = REPO_ROOT / "results" / "figures" / "boxplots_standout_genes.png"

EXCEL_CORRUPTED_GENE_SYMBOLS = ["1-Mar", "2-Mar"]
MIN_TOTAL_COUNT = 10
N_TOP_DE_GENES = 15

# Named genes of specific interest, independent of DE ranking — included
# because they're established markers in the barrier-vs-immune literature,
# even if they don't individually rank in the top DE hits.
NAMED_BARRIER_GENES = ["FLG", "LOR", "IVL", "KRT77", "LCE3D", "LCE5A", "SPRR2G"]
NAMED_IMMUNE_GENES = ["IL13", "CCL18", "TNFRSF21", "OASL", "IL4R"]

STANDOUT_GENES_FOR_BOXPLOTS = ["LCE3D", "KRT77", "CCL18", "OASL"]

CONDITION_ORDER = ["healthy", "non-lesional", "lesional"]
CONDITION_COLORS = {"healthy": "#4C9F70", "non-lesional": "#E8A33D", "lesional": "#C44E52"}


def load_data():
    counts = pd.read_csv(COUNTS_PATH, index_col=0).T
    meta = pd.read_csv(METADATA_PATH, index_col="sample_id")
    meta = meta.loc[counts.index]
    present = [g for g in EXCEL_CORRUPTED_GENE_SYMBOLS if g in counts.columns]
    if present:
        counts = counts.drop(columns=present)
    keep = counts.sum(axis=0) >= MIN_TOTAL_COUNT
    return counts.loc[:, keep], meta


def compute_vst(counts: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition")
    dds.fit_size_factors()
    dds.vst()
    return pd.DataFrame(dds.layers["vst_counts"], index=meta.index, columns=counts.columns)


def build_gene_panel(vst: pd.DataFrame) -> list:
    de = pd.read_csv(DE_PATH, index_col=0)
    top_de = de.sort_values("padj").head(N_TOP_DE_GENES).index.tolist()
    panel = list(dict.fromkeys(top_de + NAMED_BARRIER_GENES + NAMED_IMMUNE_GENES))
    present = [g for g in panel if g in vst.columns]
    missing = [g for g in panel if g not in vst.columns]
    if missing:
        print(f"WARNING: {len(missing)} panel gene(s) not found in filtered matrix: {missing}")
    return present


def plot_heatmap(vst: pd.DataFrame, meta: pd.DataFrame, genes: list):
    panel_vst = vst[genes].T  # genes x samples
    # z-score each gene across samples so the heatmap reflects relative
    # expression pattern, not absolute VST scale (which differs gene to gene)
    panel_z = panel_vst.sub(panel_vst.mean(axis=1), axis=0).div(panel_vst.std(axis=1), axis=0)

    sample_order = meta.sort_values(
        by="condition", key=lambda s: s.map({c: i for i, c in enumerate(CONDITION_ORDER)})
    ).index
    panel_z = panel_z[sample_order]

    col_colors = meta.loc[sample_order, "condition"].map(CONDITION_COLORS)

    g = sns.clustermap(
        panel_z, col_cluster=False, row_cluster=True,
        col_colors=col_colors, cmap="vlag", center=0,
        figsize=(12, 8), xticklabels=False,
        cbar_kws={"label": "z-score (VST expression)"},
    )
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xlabel(f"Samples (n={len(sample_order)}, grouped by condition)")

    handles = [plt.Rectangle((0, 0), 1, 1, color=CONDITION_COLORS[c]) for c in CONDITION_ORDER]
    g.ax_heatmap.legend(handles, CONDITION_ORDER, title="Condition",
                         bbox_to_anchor=(1.15, 1), loc="upper left", frameon=False)

    OUT_HEATMAP.parent.mkdir(parents=True, exist_ok=True)
    g.savefig(OUT_HEATMAP, dpi=150)
    print(f"Saved heatmap -> {OUT_HEATMAP}")
    plt.close()


def plot_boxplots(vst: pd.DataFrame, meta: pd.DataFrame, genes: list):
    genes = [g for g in genes if g in vst.columns]
    fig, axes = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4))
    if len(genes) == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes):
        df = pd.DataFrame({
            "expression": vst[gene],
            "condition": meta["condition"],
        })
        sns.boxplot(
            data=df, x="condition", y="expression", order=CONDITION_ORDER,
            palette=CONDITION_COLORS, ax=ax,
        )
        sns.stripplot(
            data=df, x="condition", y="expression", order=CONDITION_ORDER,
            color="black", size=3, alpha=0.5, ax=ax,
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.set_ylabel("VST expression" if ax is axes[0] else "")

    plt.tight_layout()
    OUT_BOXPLOTS.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUT_BOXPLOTS, dpi=150)
    print(f"Saved boxplots -> {OUT_BOXPLOTS}")
    plt.close()


def main():
    counts, meta = load_data()
    vst = compute_vst(counts, meta)

    genes = build_gene_panel(vst)
    print(f"Curated panel: {len(genes)} genes")

    plot_heatmap(vst, meta, genes)
    plot_boxplots(vst, meta, STANDOUT_GENES_FOR_BOXPLOTS)


if __name__ == "__main__":
    main()