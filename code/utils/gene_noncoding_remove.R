#!/usr/bin/env Rscript

# Auto-install data.table package if not available
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}
library(data.table)

# ================= Configuration Area =================
input_file <- "/shares/rheumatologie.usz/zly/proj1_2_lda_evaluation/snp_gene_projection/snp_gene_map_merged.txt"
output_file <- "/shares/rheumatologie.usz/zly/proj1_2_lda_evaluation/snp_gene_projection/snp_gene_map_merged_coding_only.txt"
hgnc_local_cache <- "hgnc_complete_set.txt" # Local cache filename

# ================= 1. Connect and fetch external authoritative database =================
cat("Connecting to HGNC external database to download latest gene annotations...\n")
# [Fixed]: Use the latest official HGNC Google Cloud Storage download URL
url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

# Use local cache if available to speed up subsequent runs
if (!file.exists(hgnc_local_cache)) {
  download.file(url, destfile = hgnc_local_cache, method = "auto", quiet = FALSE)
  cat("Database download complete!\n")
} else {
  cat("Found local HGNC database cache, loading directly...\n")
}

# Read database (quote = "" prevents parsing errors from irregular quotes within HGNC)
cat("Building protein-coding gene dictionary...\n")
hgnc_df <- fread(hgnc_local_cache, sep = "\t", quote = "")

# Filter genes explicitly classified as "protein-coding gene"
pc_genes <- hgnc_df[locus_group == 'protein-coding gene']

# Extract official symbol and remove empty values
symbols <- pc_genes$symbol[!is.na(pc_genes$symbol) & pc_genes$symbol != ""]

# Extract alias symbols and split by "|"
aliases <- pc_genes$alias_symbol[!is.na(pc_genes$alias_symbol) & pc_genes$alias_symbol != ""]
aliases_split <- unlist(strsplit(aliases, "\\|"))

# Extract previous symbols and split by "|"
prevs <- pc_genes$prev_symbol[!is.na(pc_genes$prev_symbol) & pc_genes$prev_symbol != ""]
prevs_split <- unlist(strsplit(prevs, "\\|"))

# [Core]: Merge all names into a unique set (equivalent to Python's Set)
valid_coding_genes <- unique(trimws(c(symbols, aliases_split, prevs_split)))
valid_coding_genes <- valid_coding_genes[valid_coding_genes != ""]

cat(sprintf("Successfully built! Indexed %d valid protein-coding gene names (including aliases).\n", length(valid_coding_genes)))

# ================= 2. Clean your data =================
cat("\nCleaning your SNP-Gene mapping file...\n")

# Use data.table for fast file reading (auto-split by space/tab into V1, V2 columns)
# fill=TRUE prevents errors when some SNPs have no gene following
dt <- fread(input_file, header = FALSE, col.names = c("SNP", "GENES"), fill = TRUE)
count_total <- nrow(dt)

cat("Parsing and expanding gene lists...\n")
# Split comma-separated gene names into multiple rows (e.g. one row with A,B becomes two: SNP A and SNP B)
dt_long <- dt[, .(GENE = unlist(strsplit(GENES, ","))), by = SNP]

# Trim whitespace around gene names
dt_long[, GENE := trimws(GENE)]

cat("Applying HGNC protein-coding gene filter...\n")
# [Core Logic]: Vectorized fast filtering, keep only genes present in the HGNC dictionary
dt_filtered <- dt_long[GENE %in% valid_coding_genes]

count_kept <- nrow(dt_filtered)

# Save the cleaned results
cat("Saving results...\n")
fwrite(dt_filtered, output_file, sep = " ", col.names = FALSE, quote = FALSE)

cat("========================================\n")
cat("Cleaning complete!\n")
cat(sprintf("Total rows in original file: %d\n", count_total))
cat(sprintf("Protein-coding gene mappings retained: %d\n", count_kept))
cat(sprintf("Results saved to: %s\n", output_file))