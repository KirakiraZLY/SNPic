# downstream_gaussian.R
run_downstream_gaussian <- function(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, seed, keep_all_traits=FALSE) {
  cat("\n############################################\n")
  cat("Starting Downstream Analysis (Gaussian LDA)\n")
  cat("############################################\n")

  k_topic <- best_k
  prefix_gaussian <- paste0(prefix, "_downstream_gaussian_k", k_topic)
  prefix_pathway  <- paste0(dirname(prefix), "/pathway_results/gaussian_k", k_topic)
  prefix_tissue   <- paste0(dirname(prefix), "/tissue_analysis/gaussian_k", k_topic) 

  if (!dir.exists(dirname(prefix_pathway))) dir.create(dirname(prefix_pathway), recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(prefix_tissue)) dir.create(prefix_tissue, recursive = TRUE, showWarnings = FALSE) 

  meta_df$Disease <- trimws(meta_df$Disease)
  passed_ids <- if(keep_all_traits) names(final_res$confidence_score) else names(final_res$confidence_score)[final_res$confidence_score >= best_thresh] 
  stable_list_df <- meta_df[meta_df$Disease %in% trimws(passed_ids), ]

  rownames(input_mat) <- trimws(rownames(input_mat))
  matrix_result <- input_mat[intersect(rownames(input_mat), stable_list_df$Disease), , drop = FALSE]
  matrix_result <- matrix_result[, apply(matrix_result, 2, function(x) length(unique(x)) > 1), drop = FALSE]

  # Gaussian Execution (Mixed Membership)
  dtm_scaled <- as.data.frame(t(scale(t(matrix_result))))
  dtm_scaled <- dtm_scaled - min(dtm_scaled) 

  #set.seed(seed)
  #results_umap_pca_gaussian <- plot_umap_pca_from_matrix(shared_snp_matrix = dtm_scaled, prefix = prefix_gaussian, master_map = master_map)
  set.seed(seed)
  out_gaussian <- run_mixed_membership(matrix_result, k = k_topic, master_map = master_map)
  gaussian_topic_matrix_std <- out_gaussian$topic_matrix2

  # Dimensionality reduction plots
  create_dim_plot <- function(df, x_col, y_col, title, x_lab, y_lab) {
    if(!"label2" %in% colnames(df)) df$label2 <- "Unknown"
    shape_map <- setNames(rep(c(16, 17, 15, 18, 3, 4, 8, 1, 2, 5), length.out=length(unique(df$label2))), unique(df$label2))
    p <- ggplot(df, aes_string(x=x_col, y=y_col, color="label", shape="label2", label="Disease")) +
      geom_point(size=3, alpha=0.8) + geom_text_repel(size=3.5, max.overlaps=20, bg.color="white", bg.r=0.1) +
      theme_minimal() + labs(title=title, x=x_lab, y=y_lab, color="Category", shape="Source") +
      scale_shape_manual(values=shape_map)
    return(p)
  }

  #ggsave(paste0(prefix_gaussian, "_pca.png"), create_dim_plot(results_umap_pca_gaussian$pca_df, "PC1", "PC2", "PCA (Gaussian Model)", "PC1", "PC2"), width=14, height=8, dpi=300)
  #ggsave(paste0(prefix_gaussian, "_umap.png"), create_dim_plot(results_umap_pca_gaussian$umap_df, "UMAP1", "UMAP2", "UMAP (Gaussian Model)", "UMAP1", "UMAP2"), width=14, height=8, dpi=300)

  # Similarity Analysis
  if(!"label2" %in% colnames(gaussian_topic_matrix_std)) gaussian_topic_matrix_std$label2 <- gaussian_topic_matrix_std$label
  sim_res_gaussian <- analyze_disease_similarity(gaussian_topic_matrix_std, model_type = "Gaussian", mode = "separate", edge_threshold = 0.70, seed = seed)
  if(!is.null(sim_res_gaussian$correlation_plots$network)) ggsave(paste0(prefix_gaussian, "_similarity_network.png"), sim_res_gaussian$correlation_plots$network, width=14, height=10)
  #if(!is.null(sim_res_gaussian$correlation_plots$heatmap)) ggsave(paste0(prefix_gaussian, "_similarity_heatmap.png"), sim_res_gaussian$correlation_plots$heatmap, width=14, height=10)

  # ==========================================
  # 【Core Modification Area 1】：Topic Distributions plotting function
  # ==========================================
  plot_topic_dist <- function(topic_mat, title_str, out_file) {
    non_topic_cols <- c("Disease", "disease_clean", "label_map", "label", "label2", "Disease_file", "Disease_Unique", "MaxTopic")
    topic_cols <- setdiff(colnames(topic_mat), non_topic_cols)
    
    cols_to_keep <- c(topic_cols, "Disease", if ("label2" %in% colnames(topic_mat)) "label2" else NULL)
    df_topics <- topic_mat[, cols_to_keep]
    df_topics$MaxTopic <- factor(apply(df_topics[, topic_cols], 1, which.max), levels = 1:length(topic_cols))
    
    # Handle Display_Name for both Gene-as-word and SS-as-word
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
      geom_bar(stat="identity", position="stack", width = 0.85) + # Slightly thinner bars for better visual breathing room
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
        plot.margin = margin(20, 30, 20, 20) # Increase margins to prevent label clipping
      )
      
    # 【Core Fix】：Increase width from 14 to 22 to fix horizontal compression
    ggsave(out_file, p, width=22, height=max(10, nrow(df_topics) * 0.5), limitsize = FALSE)
  }  
  plot_topic_dist(gaussian_topic_matrix_std, "Gaussian Topic Distribution", paste0(prefix_gaussian, "_topic_dist_sorted.png"))

  # ==========================================
  # 【Core Modification Area 2】：New function to plot Top Genes Per Topic
  # ==========================================
  plot_top_genes_gaussian <- function(gene_terms_mat, title_str, out_file, top_n = 10) {
    # 1. Extract Top N genes from matrix, convert to long format
    df_long <- data.frame()
    
    for (i in 1:ncol(gene_terms_mat)) {
      topic_name <- colnames(gene_terms_mat)[i]
      if (!grepl("Topic", topic_name, ignore.case=TRUE)) {
        topic_name <- paste0("Topic ", i)
      }
      
      topic_weights <- gene_terms_mat[, i]
      # Sort and extract indices with highest weights
      top_idx <- order(topic_weights, decreasing = TRUE)[1:top_n]
      
      temp_df <- data.frame(
        Topic = topic_name,
        Gene = rownames(gene_terms_mat)[top_idx],
        Weight = topic_weights[top_idx],
        stringsAsFactors = FALSE
      )
      df_long <- rbind(df_long, temp_df)
    }
    
    # 2. Normalize and generate unique IDs for plotting
    plot_terms <- df_long %>%
      group_by(Topic) %>%
      mutate(Normalized_Weight = Weight / max(Weight)) %>%
      ungroup() %>%
      mutate(unique_term_id = paste(Topic, Gene, sep = "___"))
    
    # 3. Dynamically set layout and canvas size (max 4 columns)
    n_topics <- length(unique(plot_terms$Topic))
    n_cols <- min(4, n_topics)
    n_rows <- ceiling(n_topics / n_cols)
    
    # 4. Plotting
    p_top <- ggplot(plot_terms, aes(x = reorder(unique_term_id, Normalized_Weight), y = Normalized_Weight, fill = Topic)) +
      geom_bar(stat = "identity") +
      facet_wrap(~Topic, scales = "free_y", ncol = n_cols) +
      coord_flip() +
      theme_classic() +
      scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = c(0.01, 0)) +
      scale_x_discrete(labels = function(x) sub(".*___", "", x)) +
      theme(
        text = element_text(size = 18),
        axis.text.y = element_text(size = 16, color = "black", face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        strip.text = element_text(size = 18, face = "bold"),
        legend.position = "none",
        panel.spacing.x = unit(3.0, "lines"),
        panel.spacing.y = unit(2.0, "lines")
      ) +
      labs(title = title_str, x = "", y = "Relative Importance")
    
    ggsave(out_file, p_top, width = max(14, n_cols * 4.5), height = max(8, n_rows * 4), dpi = 300, limitsize = FALSE)
  }
  
  # Call the function above, extract top 10 genes
  plot_top_genes_gaussian(out_gaussian$gene_terms, 
                          title_str = "Top 10 Genes per Topic (Gaussian Model)", 
                          out_file = paste0(prefix_gaussian, "_top10_genes.png"), 
                          top_n = 10)

  # Downstream Analyses
  tryCatch({
    run_pathway_analysis(normal_top_genes = NULL, gaussian_top_genes = out_gaussian$gene_terms,
                         normal_topic_distribution = NULL, gaussian_topic_distribution = out_gaussian$topic_matrix,
                         prefix_pathway = prefix_pathway, top_genes = 50, top_terms = 5)
    cat("Pathway analysis completed.\n")
  }, error = function(e) {cat("Pathway analysis skipped/failed:", e$message, "\n")})

  tryCatch({
    run_tissue_expression_analysis(top_genes_df = out_gaussian$gene_terms, output_dir = prefix_tissue) 
    cat("GTEx Tissue analysis completed.\n")
  }, error = function(e) {cat("GTEx analysis skipped/failed:", e$message, "\n")})

  for (f in list.files(output_dir, pattern = "SSword_", full.names = TRUE)) file.rename(f, gsub("SSword", "", f))
  cat("Gaussian LDA Analysis Complete.\n")
}