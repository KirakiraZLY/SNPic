#!/usr/bin/env Rscript

# 如果没有 data.table 包则自动安装
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}
library(data.table)

# ================= 配置区域 =================
input_file <- "/shares/rheumatologie.usz/zly/proj1_2_lda_evaluation/snp_gene_projection/snp_gene_map_merged.txt"
output_file <- "/shares/rheumatologie.usz/zly/proj1_2_lda_evaluation/snp_gene_projection/snp_gene_map_merged_coding_only.txt"
hgnc_local_cache <- "hgnc_complete_set.txt" # 本地缓存文件名

# ================= 1. 连接并获取外部权威数据库 =================
cat("正在连接 HGNC 外部数据库下载最新基因注释...\n")
# 【已修复】：使用 HGNC 最新的 Google Cloud Storage 官方下载地址
url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

# 为了加快后续运行速度，如果本地已有缓存则直接读取
if (!file.exists(hgnc_local_cache)) {
  download.file(url, destfile = hgnc_local_cache, method = "auto", quiet = FALSE)
  cat("数据库下载完成！\n")
} else {
  cat("发现本地 HGNC 数据库缓存，直接加载...\n")
}

# 读取数据库 (quote = "" 防止 HGNC 内部包含的不规则引号导致解析错位)
cat("正在构建蛋白编码基因 (Protein-Coding Genes) 词典...\n")
hgnc_df <- fread(hgnc_local_cache, sep = "\t", quote = "")

# 筛选出明确分类为 "protein-coding gene" 的基因
pc_genes <- hgnc_df[locus_group == 'protein-coding gene']

# 提取正式 symbol，并剔除空值
symbols <- pc_genes$symbol[!is.na(pc_genes$symbol) & pc_genes$symbol != ""]

# 提取别名 (alias_symbol) 并按照 "|" 分割
aliases <- pc_genes$alias_symbol[!is.na(pc_genes$alias_symbol) & pc_genes$alias_symbol != ""]
aliases_split <- unlist(strsplit(aliases, "\\|"))

# 提取曾用名 (prev_symbol) 并按照 "|" 分割
prevs <- pc_genes$prev_symbol[!is.na(pc_genes$prev_symbol) & pc_genes$prev_symbol != ""]
prevs_split <- unlist(strsplit(prevs, "\\|"))

# 【核心】：使用 unique 将所有名字合并进哈希向量 (相当于 Python 的 Set)
valid_coding_genes <- unique(trimws(c(symbols, aliases_split, prevs_split)))
valid_coding_genes <- valid_coding_genes[valid_coding_genes != ""]

cat(sprintf("成功构建！共录入 %d 个有效的蛋白编码基因名称（含别名）。\n", length(valid_coding_genes)))

# ================= 2. 清洗您的数据 =================
cat("\n开始清洗您的 SNP-Gene 映射文件...\n")

# 使用 data.table 高速读取输入文件 (按空格/制表符自动拆分为 V1, V2 列)
# fill=TRUE 防止有些 SNP 后面没有基因导致报错
dt <- fread(input_file, header = FALSE, col.names = c("SNP", "GENES"), fill = TRUE)
count_total <- nrow(dt)

cat("正在解析并展开基因列表...\n")
# 将以逗号分隔的基因名拆分成多行 (例如一行 SNP 对应 A,B 会变成两行 SNP A 和 SNP B)
dt_long <- dt[, .(GENE = unlist(strsplit(GENES, ","))), by = SNP]

# 清除基因名称前后的空白字符
dt_long[, GENE := trimws(GENE)]

cat("正在应用 HGNC 蛋白编码基因过滤器...\n")
# 【核心逻辑】：向量化极速过滤，保留存在于 HGNC 字典中的基因
dt_filtered <- dt_long[GENE %in% valid_coding_genes]

count_kept <- nrow(dt_filtered)

# 保存清洗后的结果
cat("正在保存结果...\n")
fwrite(dt_filtered, output_file, sep = " ", col.names = FALSE, quote = FALSE)

cat("========================================\n")
cat("清洗完成！\n")
cat(sprintf("处理的总行数 (原始文件行数): %d\n", count_total))
cat(sprintf("保留的蛋白编码基因映射对数: %d\n", count_kept))
cat(sprintf("结果已保存至: %s\n", output_file))