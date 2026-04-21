library(maptpx)
library(reshape2)
library(ggplot2)
library(dplyr)
library(umap)
library(ggrepel)
library(dendextend)
library(topicmodels)
library(tm)
library(pheatmap)
set.seed(123)
calculate_shared_snp_matrix <- function(folder1, folder2 = folder1, plot_heatmap = TRUE, plot_dendrogram = TRUE, n_cluster = 4, prefix = NULL) { # folder2==folder1 as default
  
  files1 <- list.files(folder1, full.names = TRUE)
  files2 <- list.files(folder2, full.names = TRUE)
  
  result <- matrix(NA, nrow = length(files1), ncol = length(files2))
  rownames(result) <- basename(files1)
  colnames(result) <- basename(files2)
  
  for (i in seq_along(files1)) {
    d1 <- read.table(files1[i], header = FALSE)
    values1 <- d1[[1]]
    print(i)
    for (j in seq_along(files2)) {
      d2 <- read.table(files2[j], header = FALSE)
      values2 <- d2[[1]]
      count <- sum(values1 %in% values2)
      result[i, j] <- count
    }
  }
  
  # Replace diagonal with maximum of other values (excluding self-comparison)
  if (folder1 == folder2) {
    for (i in 1:nrow(result)) {
      result[i, i] <- max(result[i, -i], na.rm = TRUE)
    }
  }
  
  # Remove all-zero rows and columns
  result <- result[apply(result, 1, function(x) !all(x == 0)), ]
  result <- result[, apply(result, 2, function(x) !all(x == 0))]
  
  # Plot heatmap
  if (plot_heatmap) {
    df_melt <- reshape2::melt(as.matrix(result))
    p1 <- ggplot(df_melt, aes(Var2, Var1, fill = value)) +
      geom_tile(color = NA) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(x = "Sumstats", y = "Sumstats", fill = "Shared SNPs") +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    cat("================ Plot heatmap ================ \n")
    print(p1)
    if (!is.null(prefix)) {
      ggsave(paste0(prefix, "_heatmap.png"), p1, width = 10, height = 8)
    }
  }
  
  # Hierarchical clustering
  if (plot_dendrogram) {
    result_t <- t(result)
    dist_result <- as.dist(1 - cor(result_t, method = "pearson"))
    hc_individuals <- hclust(dist_result, method = "ward.D2")
    
    pheatmap(result, cluster_rows = hc_individuals, cluster_cols = FALSE, show_colnames = FALSE, fontsize_row = 5)
    
    plot(hc_individuals, main = "Hierarchical Clustering of Sumstats", xlab = "", sub = "", cex = 0.6, las = 2)
    rect.hclust(hc_individuals, k = n_cluster, border = 1:n_cluster)
    
    if (!is.null(prefix)) {
      png(paste0(prefix, "_dendrogram.png"), width = 800, height = 600)
      plot(hc_individuals, main = "Hierarchical Clustering of Sumstats", xlab = "", sub = "", cex = 0.6, las = 2)
      rect.hclust(hc_individuals, k = n_cluster, border = 1:n_cluster)
      dev.off()
    }
  }
  
  return(result)
}

run_gene_as_word_analysis <- function(folder, snp_gene_map_file, prefix = NULL, top_n_genes = 1000) {
  require(data.table)
  require(Matrix)
  time_start <- Sys.time()
  
  # 1. 使用 data.table 极速读取映射表
  snp_gene_map <- fread(snp_gene_map_file, header = FALSE, col.names = c("SNP", "Gene"))
  setkey(snp_gene_map, SNP) 
  
  # 2. 批量读取所有表型文件
  files <- list.files(folder, full.names = TRUE)
  cat(sprintf("Processing %d files using data.table...\n", length(files)) )
  
  dt_list <- lapply(files, function(f) {
    d <- fread(f, header = FALSE, select = 1, col.names = "SNP")
    d[, file_name := basename(f)]
    return(d)
  })
  
  dt_all <- rbindlist(dt_list)
  dt_mapped <- merge(dt_all, snp_gene_map, by = "SNP", all.x = FALSE)
  dt_mapped <- dt_mapped[!grepl("XXbac", Gene, ignore.case = TRUE)]
  dt_counts <- dt_mapped[, .N, by = .(file_name, Gene)]
  
  gene_totals <- dt_counts[, .(Total = sum(N)), by = Gene][order(-Total)]
  top_genes <- head(gene_totals$Gene, top_n_genes)
  dt_counts <- dt_counts[Gene %in% top_genes]
  
  file_factors <- factor(dt_counts$file_name)
  gene_factors <- factor(dt_counts$Gene)
  
  sparse_mat <- sparseMatrix(
    i = as.integer(file_factors),
    j = as.integer(gene_factors),
    x = dt_counts$N,
    dimnames = list(levels(file_factors), levels(gene_factors))
  )
  
  result <- as.data.frame(as.matrix(sparse_mat))
  
  time_end <- Sys.time()
  cat("Gene-as-word Processing Time:", difftime(time_end, time_start, units = "secs"), "secs\n\n")
  
  return(result)
}

# 统一且纯净的 resolve_labels 函数
resolve_labels <- function(df, master_map = NULL) {
  if (is.null(master_map)) {
    df$label <- "Other"
    df$label2 <- "Unknown"
    return(df)
  }
  
  # 清理疾病名称用于匹配
  if (!"disease_clean" %in% colnames(df)) {
    df$disease_clean <- gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$", "", df$Disease)
  }
  
  # 防止重复匹配导致数据膨胀
  master_sub <- master_map[!duplicated(master_map$clean_filename), c("clean_filename", "trait", "label1", "label2")]
  
  # 合并
  df_merged <- merge(df, master_sub, by.x = "disease_clean", by.y = "clean_filename", all.x = TRUE)
  
  # 覆盖原列：有 trait 用 trait，有 label1 用 label1，否则给默认值
  df_merged$Disease <- ifelse(!is.na(df_merged$trait) & df_merged$trait != "", df_merged$trait, df_merged$Disease)
  df_merged$label <- ifelse(!is.na(df_merged$label1) & df_merged$label1 != "", df_merged$label1, "Other")
  df_merged$label2 <- ifelse(!is.na(df_merged$label2) & df_merged$label2 != "", df_merged$label2, "Unknown")
  
  # 剔除辅助列
  df_merged$disease_clean <- NULL
  df_merged$trait <- NULL
  df_merged$label1 <- NULL
  
  return(df_merged)
}

plot_umap_pca_from_matrix <- function(shared_snp_matrix, prefix=NULL, master_map = NULL, umap_n_neighbors = 5, umap_min_dist = 0.1, umap_seed = 123) {
  # shared_snp_matrix <- matrix_result
  # k = 9
  # map_disease_individual_detail = function(df) sub("^pheno_(group[0-9]+_[^\\.]+)\\..*", "\\1", df$Disease)
  # map_detail_func = map_disease_individual_detail
  # map_detail_func = function(df) sub("^pheno_(group[0-9]+_[^\\.]+)\\..*", "\\1", df$Disease)
  
  stopifnot(requireNamespace("umap"), requireNamespace("ggplot2"), requireNamespace("ggrepel"))
  
  # ============ UMAP ============
  message("Running UMAP...")
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- umap_n_neighbors
  umap_config$min_dist <- umap_min_dist
  umap_config$random_state <- umap_seed
  umap_result <- umap::umap(shared_snp_matrix, config = umap_config)
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Disease <- rownames(shared_snp_matrix)
  # label
  umap_df <- resolve_labels(umap_df, master_map)
  
  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label, label = Disease)) +
    geom_point(size = 2.5) +
    geom_text(size = 2.5, vjust = 1.5, hjust = 0.5) +
    theme_classic() +
    ggtitle("SS as word, UMAP, by label") +
    theme(legend.position = "right")
  if (!is.na(prefix)){
    ggsave(paste0(prefix, "_SSword_umap_only.png"), p_umap, width = 10, height = 8)
  }
  print(p_umap)
  
  no_label_diseases <- umap_df %>%
    filter(label == "NoLabel") %>%
    pull(Disease)
  if (length(no_label_diseases) > 0) {
    message("Diseases with NoLabel:")
    print(no_label_diseases)
  }
  
  # ============ PCA ============
  message("Running PCA...")
  pca_matrix <- prcomp(shared_snp_matrix, center = TRUE, scale. = TRUE)
  pca_df_only <- as.data.frame(pca_matrix$x[, 1:2])
  colnames(pca_df_only) <- c("PC1", "PC2")
  pca_df_only$Disease <- rownames(shared_snp_matrix)
  # label
  pca_df_only <- resolve_labels(pca_df_only, master_map)
  
  p_pca <- ggplot(pca_df_only, aes(PC1, PC2, color = label)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Disease), size = 3, max.overlaps = Inf, force = 10) +
    theme_classic() +
    ggtitle("SS as word, PCA only, by label")
  if (!is.na(prefix)){
    ggsave(paste0(prefix, "_SSword_pca_only.png"), p_pca, width = 10, height = 8)
  }
  print(p_pca)
  
  # return plots and intermediate data frames
  return(list(
    umap_plot = p_umap,
    pca_plot = p_pca,
    umap_df = umap_df,
    pca_df = pca_df_only
  ))
}


run_lda_analysis <- function(result_matrix, prefix = NULL, master_map = NULL, k = 5) {
  
  # Data Cleaning
  result <- result_matrix[apply(result_matrix, 1, function(x) !all(x == 0)), ]
  result <- result[, apply(result, 2, function(x) !all(x == 0))]
  dtm <- as.DocumentTermMatrix(as.matrix(result), weighting = weightTf)
  
  cat(" Run LDA \n")
  time_start <- Sys.time()
  
  # Run LDA
  lda_model <- LDA(
    dtm, k,
    method = "Gibbs",
    control = list(seed = 1234, alpha = 0.01, delta = 0.01, burnin = 1000, iter = 3000, thin = 100) # best = FALSE
  )
  cat(" LDA process done \n")
  
  topic_matrix <- topicmodels::posterior(lda_model)$topics
  terms_matrix <- topicmodels::posterior(lda_model)$terms
  
  ## ----------------------------------------------------------------
  ## Topic distribution: unsorted
  ## ----------------------------------------------------------------
  df_topics <- as.data.frame(topic_matrix)
  df_topics$Disease <- rownames(df_topics)
  df_melted <- reshape2::melt(df_topics, id.vars = "Disease", variable.name = "Topic", value.name = "Probability")
  df_melted$Disease <- gsub("^finngen_R12_|\\.list$", "", df_melted$Disease)
  
  p_topics <- ggplot(df_melted, aes(x = Probability, y = Disease, fill = Topic)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    labs(title = "Disease distribution across LDA topics", 
         x = "Topic Probability", 
         y = "Disease") +
    theme(axis.text.y = element_text(size = 12), 
          legend.position = "right")
  
  if (!is.null(prefix)) {
    ggsave(paste0(prefix, "_SSword_snpic_topic_distribution.png"), p_topics, width = 10, height = 12)
  } else {
    print(p_topics)
  }
  
  ## ----------------------------------------------------------------
  ## Topic distribution: sorted by max topic
  ## ----------------------------------------------------------------
  df_topics$MaxTopic <- apply(df_topics[, 1:k], 1, which.max)
  df_topics$MaxTopic <- factor(df_topics$MaxTopic, levels = 1:k, labels = paste0("Topic", 1:k))
  
  df_topics <- df_topics %>%
    arrange(MaxTopic, Disease) %>%
    mutate(Disease = factor(Disease, levels = Disease))
  
  df_melted <- reshape2::melt(df_topics, id.vars = c("Disease", "MaxTopic"), variable.name = "Topic", value.name = "Probability")
  
  p_sorted <- ggplot(df_melted, aes(x = Probability, y = Disease, fill = Topic)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_classic() +
    labs(title = "LDA topics (sorted by max topic)", 
         x = "Topic Probability", 
         y = "Disease") +
    theme(axis.text.y = element_text(size = 12),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(10, 10, 10, 10),
          legend.position = "right")
  
  print(p_sorted)
  if (!is.null(prefix)) {
    ggsave(paste0(prefix, "_SSword_snpic_topic_distribution_sorted.png"), p_sorted, width = 10, height = 12)
  }
  
  ## ----------------------------------------------------------------
  ## Top words (Modified Section with Normalization)
  ## ----------------------------------------------------------------
  df_terms <- as.data.frame(t(terms_matrix))
  df_terms$Comorbidity <- rownames(df_terms)
  df_terms_melted <- reshape2::melt(df_terms, id.vars = "Comorbidity", variable.name = "Topic", value.name = "Weight")
  
  # 1. Get Top 50 for the return object
  top_terms <- df_terms_melted %>% 
    group_by(Topic) %>% 
    arrange(desc(Weight)) %>% 
    slice_head(n = 50) %>%
    ungroup()
  
  top_terms$Comorbidity <- gsub("^finngen_R12_|\\.list$", "", top_terms$Comorbidity)
  
  # 2. Prepare Data for Plotting: Extract Top 10 AND Normalize
  plot_terms <- top_terms %>%
    group_by(Topic) %>%
    arrange(desc(Weight)) %>% 
    slice_head(n = 10) %>% 
    # --- [关键修改 1]：归一化权重，使每个 Topic 的最大值都为 1 ---
    mutate(Normalized_Weight = Weight / max(Weight)) %>%
    # --------------------------------------------------------
  ungroup() %>%
    mutate(unique_term_id = paste(Topic, Comorbidity, sep = "___"))
  
  # 3. Plotting
  # 注意：这里 y 轴变成了 Normalized_Weight
  p_top_words <- ggplot(plot_terms, aes(x = reorder(unique_term_id, Normalized_Weight), y = Normalized_Weight, fill = Topic)) +
    geom_bar(stat = "identity") +
    # scales="free_y" 允许不同分面有不同的词（y轴），但 x轴（数值）现在都是 0-1
    facet_wrap(~Topic, scales = "free_y") + 
    coord_flip() +
    theme_classic() +
    
    # --- [关键修改 2]：控制横轴（原Y轴）刻度 ---
    scale_y_continuous(
      limits = c(0, 1.05), # 强制范围 0 到 1 (多给0.05防止顶格)
      breaks = c(1),       # 只显示最大值 1
      expand = c(0, 0)     # 减少边缘留白
    ) +
    # ------------------------------------------
  
  # 去除 label 前缀
  scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
    
    theme(
      text = element_text(size = 24),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 24),
      plot.title = element_text(size = 24),
      strip.text = element_text(size = 24)
    ) +
    labs(title = "Top 10 Words per Topic (Normalized)", x = "Word", y = "Relative Importance")
  
  if (!is.null(prefix)) {
    ggsave(paste0(prefix, "_snpic_top10_words_normalized.png"), p_top_words, width = 18, height = 16)
  }
  print(p_top_words)
  
  ## ----------------------------------------------------------------
  ## UMAP
  ## ----------------------------------------------------------------
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- 5
  umap_config$random_state <- 123
  
  umap_res <- umap(topic_matrix, config = umap_config)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Disease <- rownames(topic_matrix)
  umap_df$MaxTopic <- df_topics$MaxTopic
  
  # Label
  umap_df <- resolve_labels(umap_df, master_map)
  
  # Get colors
  custom_colors <- get_network_colors()
  unique_labels_umap <- unique(umap_df$label)
  n_labels_umap <- length(unique_labels_umap)
  
  p_snpic_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = label)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Disease), size = 3, max.overlaps = Inf) +
    theme_classic() +
    ggtitle("SS as word, SNPic-SS + UMAP Clustering, by label") +
    scale_color_manual(
      values = if(n_labels_umap <= length(custom_colors)) {
        setNames(custom_colors[1:n_labels_umap], unique_labels_umap)
      } else {
        NULL
      }
    ) +
    theme(legend.position = "right")
  
  if (!is.null(prefix)) {
    ggsave(paste0(prefix, "_SSword_snpic_umap.png"), p_snpic_umap, width = 10, height = 8)
  }
  print(p_snpic_umap)
  
  ## ----------------------------------------------------------------
  ## PCA
  ## ----------------------------------------------------------------
  pca_res <- prcomp(topic_matrix, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$Disease <- rownames(topic_matrix)
  pca_df$MaxTopic <- df_topics$MaxTopic
  
  # Label
  pca_df <- resolve_labels(pca_df, master_map)
  
  unique_labels_pca <- unique(pca_df$label)
  n_labels_pca <- length(unique_labels_pca)
  
  p_snpic_pca <- ggplot(pca_df, aes(PC1, PC2, color = label)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Disease), size = 3, max.overlaps = Inf) +
    theme_classic() +
    ggtitle("SS as word, SNPic-SS + PCA Clustering, by label") +
    scale_color_manual(
      values = if(n_labels_pca <= length(custom_colors)) {
        setNames(custom_colors[1:n_labels_pca], unique_labels_pca)
      } else {
        NULL
      }
    ) +
    theme(legend.position = "right")
  
  if (!is.null(prefix)) {
    ggsave(paste0(prefix, "_SSword_snpic_pca.png"), p_snpic_pca, width = 10, height = 8)
  }
  print(p_snpic_pca)
  
  time_end <- Sys.time()
  cat("LDA + PCA Processing Time: ", time_end - time_start, "\n\n")
  
  no_label_diseases <- pca_df %>% filter(label == "NoLabel") %>% pull(Disease)
  print(no_label_diseases)
  
  topic_matrix2 <- as.data.frame(topic_matrix)
  topic_matrix2$Disease <- rownames(topic_matrix)
  topic_matrix2 <- resolve_labels(topic_matrix2, master_map)
  
  # Return top 50 in the list as requested, but we plotted top 10
  return(list(lda_model = lda_model, 
              topic_matrix = topic_matrix2, 
              top_word = top_terms, 
              umap_df = umap_df, 
              pca_df = pca_df))
}

# 首先，创建一个函数来计算困惑度和一致性
calculate_lda_metrics <- function(result_matrix, k) {
  # 清理输入矩阵
  result <- result_matrix[apply(result_matrix, 1, function(x) !all(x == 0)), ]
  result <- result[, apply(result, 2, function(x) !all(x == 0))]
  
  # 创建文档-词项矩阵
  dtm <- as.DocumentTermMatrix(as.matrix(result), weighting = weightTf)
  
  # 运行LDA
  lda_model <- topicmodels::LDA(
    dtm, k,
    method = "Gibbs",
    control = list(seed = 1234, alpha = 0.01, delta = 0.01, burnin = 1000, iter = 3000, thin = 100)
  )
  
  # 计算困惑度
  perp <- perplexity(lda_model, newdata = dtm)
  
  # 计算一致性 (简化版本)
  phi <- exp(lda_model@beta)  # 主题-词分布
  
  # 计算主题间相似性作为一致性代理
  topic_similarity <- cor(t(phi))
  diag(topic_similarity) <- NA  # 忽略对角线
  coherence <- mean(topic_similarity, na.rm = TRUE)
  
  return(list(perplexity = perp, coherence = coherence))
}

choose_best_n_topic_multi <- function(result_matrix, 
                                  prefix = NULL, 
                                  map_disease_individual_detail = map_disease_individual_detail,
                                  map_detail_func = map_disease_individual_detail,
                                  topic_range = 4:10) {
  
  # 评估不同主题数量
  metrics <- data.frame(k = topic_range, perplexity = NA, coherence = NA)
  
  for (k in topic_range) {
    cat("Evaluating k =", k, "\n")
    metrics_result <- calculate_lda_metrics(result_matrix, k)
    metrics[metrics$k == k, "perplexity"] <- metrics_result$perplexity
    metrics[metrics$k == k, "coherence"] <- metrics_result$coherence
  }
  
  # 选择最佳k值
  # 归一化指标并计算综合评分
  metrics$norm_perplexity <- 1 / metrics$perplexity  # 困惑度越低越好，所以取倒数
  metrics$norm_coherence <- metrics$coherence  # 一致性越高越好
  
  # 归一化到0-1范围
  metrics$norm_perplexity <- (metrics$norm_perplexity - min(metrics$norm_perplexity, na.rm = TRUE)) / 
    (max(metrics$norm_perplexity, na.rm = TRUE) - min(metrics$norm_perplexity, na.rm = TRUE))
  
  metrics$norm_coherence <- (metrics$norm_coherence - min(metrics$norm_coherence, na.rm = TRUE)) / 
    (max(metrics$norm_coherence, na.rm = TRUE) - min(metrics$norm_coherence, na.rm = TRUE))
  
  # 计算综合评分（各指标权重可调整）
  metrics$composite_score <- 0.6 * metrics$norm_perplexity + 0.4 * metrics$norm_coherence
  
  # 选择综合评分最高的k
  best_k <- metrics$k[which.max(metrics$composite_score)]
  cat("Selected best k:", best_k, "\n")
  
  # 使用最佳k运行完整分析
  result <- run_lda_analysis(
    result_matrix = result_matrix,
    prefix = prefix,
    # map_disease_individual_detail = map_disease_individual_detail,
    map_detail_func = map_detail_func,
    k = best_k
  )
  
  # 添加评估指标到结果
  result$best_k <- best_k
  result$metrics <- metrics
  
  return(result)
}