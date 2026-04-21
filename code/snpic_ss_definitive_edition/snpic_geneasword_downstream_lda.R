# downstream_lda.R
run_downstream_lda <- function(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, seed, keep_all_traits=FALSE) {
  cat("\n############################################\n")
  cat("Starting Downstream Analysis (Normal LDA)\n")
  cat("############################################\n")

  k_topic <- best_k
  prefix_downstream <- paste0(prefix, "_downstream_normal_k", k_topic)
  prefix_pathway  <- paste0(dirname(prefix), "/pathway_results/normal_k", k_topic)
  prefix_tissue   <- paste0(dirname(prefix), "/tissue_analysis/normal_k", k_topic) 

  if (!dir.exists(dirname(prefix_pathway))) dir.create(dirname(prefix_pathway), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(prefix_tissue)) dir.create(prefix_tissue, recursive = TRUE, showWarnings = FALSE) 

  meta_df$Disease <- trimws(meta_df$Disease)
  passed_ids <- if(keep_all_traits) names(final_res$confidence_score) else names(final_res$confidence_score)[final_res$confidence_score >= best_thresh] 
  stable_list_df <- meta_df[meta_df$Disease %in% trimws(passed_ids), ]

  rownames(input_mat) <- trimws(rownames(input_mat))
  matrix_result <- input_mat[intersect(rownames(input_mat), stable_list_df$Disease), , drop = FALSE]
  matrix_result <- matrix_result[, apply(matrix_result, 2, function(x) length(unique(x)) > 1), drop = FALSE]

  # Normal LDA Execution
  set.seed(seed)
  results_umap_pca <- plot_umap_pca_from_matrix(shared_snp_matrix = matrix_result, prefix = prefix_downstream, master_map = master_map)
  set.seed(seed)
  results_snpic <- run_lda_analysis(matrix_result, k = k_topic, prefix = prefix_downstream, master_map = master_map)

  snpic_topic_matrix <- results_snpic$topic_matrix
  snpic_top_word <- results_snpic$top_word

  # Dimensionality Reduction Plots
  create_dim_plot <- function(df, x_col, y_col, title, x_lab, y_lab) {
    if(!"label2" %in% colnames(df)) df$label2 <- "Unknown"
    shape_map <- setNames(rep(c(16, 17, 15, 18, 3, 4, 8, 1, 2, 5), length.out=length(unique(df$label2))), unique(df$label2))
    p <- ggplot(df, aes_string(x=x_col, y=y_col, color="label", shape="label2", label="Disease")) +
      geom_point(size=3, alpha=0.8) + geom_text_repel(size=3.5, max.overlaps=20, bg.color="white", bg.r=0.1) +
      theme_minimal() + labs(title=title, x=x_lab, y=y_lab, color="Category", shape="Source") +
      scale_shape_manual(values=shape_map)
    return(p)
  }

  ggsave(paste0(prefix_downstream, "_pca.png"), create_dim_plot(results_umap_pca$pca_df, "PC1", "PC2", "PCA (LDA Model)", "PC1", "PC2"), width=14, height=8, dpi=300)
  ggsave(paste0(prefix_downstream, "_umap.png"), create_dim_plot(results_umap_pca$umap_df, "UMAP1", "UMAP2", "UMAP (LDA Model)", "UMAP1", "UMAP2"), width=14, height=8, dpi=300)

  # Similarity Analysis
  if(!"label2" %in% colnames(snpic_topic_matrix)) snpic_topic_matrix$label2 <- snpic_topic_matrix$label
  sim_res_lda <- analyze_disease_similarity(snpic_topic_matrix, model_type = "LDA", mode = "separate", edge_threshold = 0.70, seed = seed)
  if(!is.null(sim_res_lda$correlation_plots$network)) ggsave(paste0(prefix_downstream, "_similarity_network.png"), sim_res_lda$correlation_plots$network, width=14, height=10)
  if(!is.null(sim_res_lda$correlation_plots$heatmap)) ggsave(paste0(prefix_downstream, "_similarity_heatmap.png"), sim_res_lda$correlation_plots$heatmap, width=14, height=10)

  # ==========================================
  # 【核心修改区】：Topic Distributions 画图函数
  # ==========================================
  plot_topic_dist <- function(topic_mat, title_str, out_file) {
    non_topic_cols <- c("Disease", "disease_clean", "label_map", "label", "label2", "Disease_file", "Disease_Unique", "MaxTopic")
    topic_cols <- setdiff(colnames(topic_mat), non_topic_cols)
    
    cols_to_keep <- c(topic_cols, "Disease", if ("label2" %in% colnames(topic_mat)) "label2" else NULL)
    df_topics <- topic_mat[, cols_to_keep]
    df_topics$MaxTopic <- factor(apply(df_topics[, topic_cols], 1, which.max), levels = 1:length(topic_cols))
    
    # 这里兼容 Gene-as-word 和 SS-as-word 的 Display_Name 逻辑
    if("label2" %in% colnames(df_topics)) {
        df_topics$Display_Name <- paste0(df_topics$Disease, " (", df_topics$label2, ")")
    } else {
        df_topics$Display_Name <- df_topics$Disease
    }
    
    df_topics <- df_topics %>% arrange(MaxTopic, Display_Name)
    df_topics$Disease_Unique <- factor(make.unique(as.character(df_topics$Display_Name)), levels = make.unique(as.character(df_topics$Display_Name)))
    
    df_melt <- reshape2::melt(df_topics, id.vars = intersect(colnames(df_topics), c("Disease_Unique", "Display_Name", "Disease", "MaxTopic", "label2")), measure.vars = topic_cols, variable.name = "Topic", value.name = "Probability")
    df_melt$TopicNum <- gsub("\\D", "", as.character(df_melt$Topic))
    
    p <- ggplot(df_melt, aes(x=Probability, y=Disease_Unique, fill=Topic)) +
      geom_bar(stat="identity", position="stack", width = 0.85) + # 柱子稍微变细一点点，增加透气感
      geom_text(data = subset(df_melt, Probability >= 0.05), aes(label = TopicNum), position = position_stack(vjust = 0.5), size = 6, color = "black", fontface = "bold") +
      scale_y_discrete(labels=setNames(as.character(df_topics$Display_Name), df_topics$Disease_Unique)) + 
      theme_classic() + 
      labs(title=title_str, x="Topic Probability", y="") + 
      theme(
        plot.title = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        axis.text.y = element_text(size = 18, face = "bold", color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 15)),
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 18, face = "bold"),
        plot.margin = margin(20, 30, 20, 20) # 增加四周的留白，防止标签被切掉
      )
      
    # 【核心修复】：将 width 从 14 提高到 22，彻底解决横向挤压的问题
    ggsave(out_file, p, width=22, height=max(10, nrow(df_topics) * 0.5), limitsize = FALSE)
  }
  
  plot_topic_dist(snpic_topic_matrix, "SNPic LDA Topic Distribution", paste0(prefix_downstream, "_topic_dist_sorted.png"))

  # Downstream Analyses
  tryCatch({
    run_pathway_analysis(normal_top_genes = snpic_top_word, gaussian_top_genes = NULL,
                         normal_topic_distribution = snpic_topic_matrix, gaussian_topic_distribution = NULL,
                         prefix_pathway = prefix_pathway, top_genes = 50, top_terms = 5)
    cat("Pathway analysis completed.\n")
  }, error = function(e) {cat("Pathway analysis skipped/failed:", e$message, "\n")})

  tryCatch({
    run_tissue_expression_analysis(top_genes_df = snpic_top_word, output_dir = prefix_tissue) 
    cat("GTEx Tissue analysis completed.\n")
  }, error = function(e) {cat("GTEx analysis skipped/failed:", e$message, "\n")})

  for (f in list.files(output_dir, pattern = "SSword_", full.names = TRUE)) file.rename(f, gsub("SSword", "", f))
  cat("Normal LDA Analysis Complete.\n")
}