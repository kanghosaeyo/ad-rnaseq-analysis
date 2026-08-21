# Raw data

This directory holds the raw GEO supplementary file. It is not tracked in
git because it's a large, re-downloadable third-party file.

## To obtain it

1. Go to [GSE121212 on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121212)
2. Under "Supplementary file", download `GSE121212_readcount.txt.gz`
3. Decompress and place it here as `data/raw/GSE121212_readcount.txt`

Expected shape after download: 31,364 genes x 147 samples, raw integer
counts, gene symbols as row index.