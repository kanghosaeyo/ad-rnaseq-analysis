"""
Gene set enrichment analysis on differential expression results.

Uses a ranked-list (GSEA-style) method against gene set libraries fetched
live from Enrichr, rather than a hard significance cutoff — this avoids
the boundary-sensitivity issue documented for the paired comparison (see
project notes), since ranked enrichment uses every gene's test statistic
rather than a binary significant/not-significant call.

"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import gseapy as gp

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

TABLES_DIR = REPO_ROOT / "results" / "tables"
FIGURES_DIR = REPO_ROOT / "results" / "figures"

# Enrichr library names. See https://maayanlab.cloud/Enrichr/#libraries
# for the full list if you want to add/swap libraries.
GENE_SET_LIBRARIES = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
]

MIN_SET_SIZE = 5
MAX_SET_SIZE = 500
N_PERMUTATIONS = 1000
RANDOM_SEED = 42
TOP_N_TERMS = 15 

# Keywords to flag terms relevant to the barrier-vs-immune question, so
# they're easy to spot inside a large, generic results table.
HIGHLIGHT_KEYWORDS = [
    "keratin", "cornif", "epiderm", "skin barrier", "desmosome",
    "interleukin", "cytokine", "chemokine", "th2", "t helper",
    "immune", "inflammat",
]


def load_ranked_genes(de_table_path: Path) -> pd.Series:
    res = pd.read_csv(de_table_path, index_col=0)
    res = res.dropna(subset=["stat"])
    return res["stat"].sort_values(ascending=False)


def run_enrichment(ranked_genes: pd.Series, library: str) -> pd.DataFrame:
    pre_res = gp.prerank(
        rnk=ranked_genes,
        gene_sets=library,
        min_size=MIN_SET_SIZE,
        max_size=MAX_SET_SIZE,
        permutation_num=N_PERMUTATIONS,
        seed=RANDOM_SEED,
        outdir=None,
    )
    results = pre_res.res2d.copy()
    results["NES"] = results["NES"].astype(float)
    results["FDR q-val"] = results["FDR q-val"].astype(float)
    return results.sort_values("FDR q-val")


def flag_relevant_terms(results: pd.DataFrame) -> pd.DataFrame:
    pattern = "|".join(HIGHLIGHT_KEYWORDS)
    results = results.copy()
    results["relevant_to_barrier_immune"] = results["Term"].str.contains(
        pattern, case=False, regex=True
    )
    return results


def plot_top_terms(results: pd.DataFrame, comparison_name: str, library: str, out_path: Path):
    sig = results[results["FDR q-val"] < 0.05]
    top = sig.reindex(sig["NES"].abs().sort_values(ascending=False).index).head(TOP_N_TERMS)
    if len(top) == 0:
        print(f"No significant terms (FDR < 0.05) for {library} — skipping plot")
        return
    top = top.sort_values("NES")

    colors = ["crimson" if flag else "steelblue"
              for flag in top["relevant_to_barrier_immune"]]

    plt.figure(figsize=(8, max(3, 0.35 * len(top))))
    plt.barh(top["Term"], top["NES"], color=colors)
    plt.axvline(0, color="black", linewidth=0.5)
    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.title(f"Top enriched terms: {comparison_name} ({library})")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    print(f"Saved plot -> {out_path}")
    plt.close()


def run_for_comparison(comparison_name: str, de_table_filename: str):
    de_path = TABLES_DIR / de_table_filename
    ranked_genes = load_ranked_genes(de_path)
    print(f"\n=== {comparison_name}: ranked {len(ranked_genes)} genes ===")

    for library in GENE_SET_LIBRARIES:
        print(f"\n--- Library: {library} ---")
        results = run_enrichment(ranked_genes, library)
        results = flag_relevant_terms(results)

        n_sig = (results["FDR q-val"] < 0.05).sum()
        n_relevant_sig = ((results["FDR q-val"] < 0.05) & results["relevant_to_barrier_immune"]).sum()
        print(f"Significant terms (FDR < 0.05): {n_sig} / {len(results)}")
        print(f"Of those, flagged as barrier/immune-relevant: {n_relevant_sig}")

        print("\nTop 10 relevant significant terms:")
        relevant_sig = results[(results["FDR q-val"] < 0.05) & results["relevant_to_barrier_immune"]]
        print(relevant_sig.sort_values("FDR q-val").head(10)[["Term", "NES", "FDR q-val"]])

        safe_lib_name = library.replace(" ", "_")
        table_path = TABLES_DIR / f"enrichment_{comparison_name}_{safe_lib_name}.csv"
        table_path.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(table_path, index=False)
        print(f"Saved full results table -> {table_path}")

        figure_path = FIGURES_DIR / f"enrichment_{comparison_name}_{safe_lib_name}_top.png"
        plot_top_terms(results, comparison_name, library, figure_path)


def main():
    run_for_comparison("lesional_vs_healthy", "de_lesional_vs_healthy.csv")
    run_for_comparison("nonlesional_vs_healthy", "de_nonlesional_vs_healthy.csv")
    run_for_comparison("lesional_vs_nonlesional", "de_lesional_vs_nonlesional.csv")


if __name__ == "__main__":
    main()