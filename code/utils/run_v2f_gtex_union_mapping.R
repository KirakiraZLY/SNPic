#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

# ================= 1. Configuration =================
input_file <- "/home/leyzha/data/proj1_2_lda_evaluation/snpic_ss_definitive_edition/snp_gene_map_merged.txt"
output_file <- "/home/leyzha/data/proj1_2_lda_evaluation/gtex_v2f/v2f_snp_gene_map_raw.txt"

gtex_dir <- "/home/leyzha/data/proj1_2_lda_evaluation/gtex_v2f"
lookup_file <- file.path(gtex_dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
gene_annot_file <- file.path(gtex_dir, "gencode.v26.GRCh38.genes.gtf")
eqtl_folder <- file.path(gtex_dir, "GTEx_Analysis_v8_eQTL")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# ================= 2. Read Original GRCh37 Data =================
cat("1. Reading your original GRCh37 SNP list...\n")
dt <- fread(input_file, header=FALSE, fill=TRUE, col.names=c("SNP", "Nearest_Gene"))
dt <- dt[grepl("^rs", SNP)] 
my_snps <- unique(dt$SNP)
cat(sprintf("   Loaded %d unique rsIDs to query.\n", length(my_snps)))

# ================= 3. Load GTEx GTF Dictionary =================
cat("2. Loading GTEx Gene Annotations...\n")
gtf <- fread(cmd = paste0("grep -v '^#' ", gene_annot_file), sep = "\t", header = FALSE)
gtf_genes <- gtf[V3 == "gene"]
gene_ids <- sub('.*gene_id "([^"]+)".*', '\\1', gtf_genes$V9)
gene_names <- sub('.*gene_name "([^"]+)".*', '\\1', gtf_genes$V9)
gene_map <- unique(data.table(gene_id = gene_ids, gene_name = gene_names))
rm(gtf, gtf_genes); gc()

# ================= 4. Cross-Build rsID Matching =================
cat("3. Translating IDs via GTEx Lookup Table (Bridging GRCh37 to GRCh38)...\n")
lookup <- fread(lookup_file, select = c("variant_id", "rs_id_dbSNP151_GRCh38p7"))
setnames(lookup, "rs_id_dbSNP151_GRCh38p7", "rsID")

# Match using your rsID, seamlessly ignoring coordinate system differences
lookup_filtered <- lookup[rsID %in% my_snps]
rm(lookup); gc()

cat(sprintf("   Matched %d variant_ids in GTEx using rsIDs.\n", nrow(lookup_filtered)))

# ================= 5. Multi-Tissue eQTL Scanning =================
cat("4. Scanning 54 GTEx tissues for significant eQTLs...\n")
eqtl_files <- list.files(eqtl_folder, pattern = "\\.signif_variant_gene_pairs\\.txt\\.gz$", full.names = TRUE)

all_eqtls <- list()
valid_variant_ids <- lookup_filtered$variant_id

for (i in seq_along(eqtl_files)) {
  f <- eqtl_files[i]
  tissue_name <- gsub("\\.v8\\.signif_variant_gene_pairs\\.txt\\.gz", "", basename(f))
  cat(sprintf("   [%d/%d] Scanning %s...\n", i, length(eqtl_files), tissue_name))
  
  tmp <- fread(f, select = c("variant_id", "gene_id"))
  tmp <- tmp[variant_id %in% valid_variant_ids]
  if(nrow(tmp) > 0) all_eqtls[[i]] <- tmp
}

# ================= 6. Merge & Translate =================
cat("5. Merging multi-tissue evidence...\n")
combined_eqtls <- unique(rbindlist(all_eqtls))
combined_eqtls <- merge(combined_eqtls, lookup_filtered, by = "variant_id", all.x = TRUE)
combined_eqtls <- merge(combined_eqtls, gene_map, by = "gene_id", all.x = TRUE)

final_v2f <- combined_eqtls[!is.na(rsID) & !is.na(gene_name), 
                            .(eGenes = paste(unique(gene_name), collapse = ",")), 
                            by = .(SNP = rsID)]

# ================= 7. Union Strategy =================
cat("6. Finalizing V2F mapping (Union of GRCh37 Nearest-Gene and 54-Tissue eQTLs)...\n")
final_df <- merge(dt, final_v2f, by="SNP", all.x=TRUE)

final_df[, Final_Gene := Nearest_Gene]
final_df[!is.na(eGenes) & eGenes != "", Final_Gene := paste(Nearest_Gene, eGenes, sep=",")]

cat("7. Deduplicating combined genes...\n")
final_df[, Final_Gene := sapply(strsplit(Final_Gene, ","), function(x) {
  paste(unique(trimws(x[!is.na(x) & x != ""])), collapse=",")
})]

upgraded_count <- sum(!is.na(final_df$eGenes) & final_df$eGenes != "")
cat(sprintf("\n>>> SUMMARY: Successfully expanded %d SNPs with functional eGenes across 54 Tissues! <<<\n", upgraded_count))

# ================= 8. Export =================
cat("8. Exporting to file...\n")
out_dt <- final_df[, .(SNP, Final_Gene)]
fwrite(out_dt, output_file, sep=" ", col.names=FALSE, quote=FALSE)

cat("========================================\n")
cat("DONE! Next step: Run gene_noncoding_remove.R on this raw file.\n")