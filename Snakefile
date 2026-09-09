"""
Snakemake workflow for the AD skin transcriptomics pipeline.

Run with: snakemake --cores 1
Dry-run: snakemake -n
"""

RAW = "data/raw/GSE121212_readcount.txt"

rule all:
    input:
        "results/tables/library_size_qc.csv",
        "results/tables/pca_coordinates.csv",
        "results/figures/pca_plot.png",
        "results/tables/de_lesional_vs_healthy.csv",
        "results/figures/volcano_lesional_vs_healthy.png",
        "results/figures/ma_lesional_vs_healthy.png",
        "results/tables/de_nonlesional_vs_healthy.csv",
        "results/figures/volcano_nonlesional_vs_healthy.png",
        "results/tables/de_lesional_vs_nonlesional.csv",
        "results/figures/volcano_lesional_vs_nonlesional.png",
        "results/figures/heatmap_curated_panel.png",
        "results/figures/boxplots_standout_genes.png",
        "results/.enrichment_done",

rule curate_subset:
    input:
        raw = RAW,
        script = "scripts/curate_subset.py",
    output:
        metadata = "data/processed/metadata.csv",
        counts = "data/processed/counts_subset.csv",
    shell:
        "python {input.script}"

rule qc_pca:
    input:
        metadata = "data/processed/metadata.csv",
        counts = "data/processed/counts_subset.csv",
        script = "scripts/qc_pca.py",
    output:
        libsize = "results/tables/library_size_qc.csv",
        pca_table = "results/tables/pca_coordinates.csv",
        pca_plot = "results/figures/pca_plot.png",
    shell:
        "python {input.script}"

rule differential_expression:
    input:
        metadata = "data/processed/metadata.csv",
        counts = "data/processed/counts_subset.csv",
        script = "scripts/differential_expression.py",
    output:
        table = "results/tables/de_lesional_vs_healthy.csv",
        volcano = "results/figures/volcano_lesional_vs_healthy.png",
        ma = "results/figures/ma_lesional_vs_healthy.png",
    shell:
        "python {input.script}"

rule secondary_comparison:
    input:
        metadata = "data/processed/metadata.csv",
        counts = "data/processed/counts_subset.csv",
        script = "scripts/secondary_comparison.py",
    output:
        table1 = "results/tables/de_nonlesional_vs_healthy.csv",
        volcano1 = "results/figures/volcano_nonlesional_vs_healthy.png",
        table2 = "results/tables/de_lesional_vs_nonlesional.csv",
        volcano2 = "results/figures/volcano_lesional_vs_nonlesional.png",
    shell:
        "python {input.script}"

# Enrichment produces many files (2 libraries x 3 comparisons x table+plot).
# A sentinel/flag file keeps this rule simple rather than hardcoding all
# ~12 output paths individually.
rule enrichment_analysis:
    input:
        de1 = "results/tables/de_lesional_vs_healthy.csv",
        de2 = "results/tables/de_nonlesional_vs_healthy.csv",
        de3 = "results/tables/de_lesional_vs_nonlesional.csv",
        script = "scripts/enrichment_analysis.py",
    output:
        flag = touch("results/.enrichment_done"),
    shell:
        "python {input.script}"

rule visualize_results:
    input:
        metadata = "data/processed/metadata.csv",
        counts = "data/processed/counts_subset.csv",
        de_table = "results/tables/de_lesional_vs_healthy.csv",
        script = "scripts/visualize_results.py",
    output:
        heatmap = "results/figures/heatmap_curated_panel.png",
        boxplots = "results/figures/boxplots_standout_genes.png",
    shell:
        "python {input.script}"