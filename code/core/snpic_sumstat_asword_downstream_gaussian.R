# downstream_ss_gaussian.R
run_downstream_ss_gaussian <- function(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, seed, keep_all_traits) {
  cat("\n############################################\n")
  cat("Starting Downstream Analysis (SS-as-Word: Gaussian LDA)\n")
  cat("############################################\n")

  k_topic <- best_k
  prefix_gaussian <- paste0(prefix, "_downstream_gaussian_k", k_topic)

  meta_df$Disease <- trimws(meta_df$Disease)
  passed_ids <- if (keep_all_traits) names(final_res$confidence_score) else names(final_res$confidence_score)[final_res$confidence_score >= best_thresh]
  passed_ids <- trimws(passed_ids)

  stable_list_df <- meta_df[meta_df$Disease %in% passed_ids, ]
  if (nrow(stable_list_df) < 2) stop("Error: Not enough stable diseases (<2) for downstream.")

  target_diseases <- stable_list_df$Disease
  rownames(input_mat) <- trimws(rownames(input_mat))
  common_diseases <- intersect(rownames(input_mat), target_diseases)

  matrix_result <- input_mat[common_diseases, common_diseases, drop = FALSE]
  valid_idx <- apply(matrix_result, 2, function(x) length(unique(x)) > 1)
  matrix_result <- matrix_result[valid_idx, valid_idx, drop = FALSE]

  # Run SNPic (Gaussian)
  dtm_scaled <- as.data.frame(t(scale(t(matrix_result))))
  dtm_scaled <- dtm_scaled - min(dtm_scaled, na.rm = TRUE)

  # set.seed(seed)
  # results_umap_pca_gaussian <- plot_umap_pca_from_matrix(dtm_scaled, prefix_gaussian, master_map)
  set.seed(seed)
  out_gaussian <- run_mixed_membership(matrix_result, k_topic, master_map = master_map)
  gaussian_topic_matrix_std <- out_gaussian$topic_matrix2

  # Similarity Analysis
  sim_res_gaussian <- analyze_disease_similarity(gaussian_topic_matrix_std, model_type = "Gaussian", mode = "separate", edge_threshold = 0.70, seed = seed)
  if(!is.null(sim_res_gaussian$correlation_plots$network)) ggsave(paste0(prefix_gaussian, "_similarity_network_gaussian.png"), sim_res_gaussian$correlation_plots$network, width=14, height=10)

  # ==========================================
  # 1. Topic Distributions plot function
  # ==========================================
  plot_topic_dist <- function(topic_mat, title_str, out_file) {
    non_topic_cols <- c("Disease", "disease_clean", "label_map", "label", "label2", "Disease_file", "Disease_Unique", "MaxTopic")
    topic_cols <- setdiff(colnames(topic_mat), non_topic_cols)
    
    cols_to_keep <- c(topic_cols, "Disease", if ("label2" %in% colnames(topic_mat)) "label2" else NULL)
    df_topics <- topic_mat[, cols_to_keep]
    df_topics$MaxTopic <- factor(apply(df_topics[, topic_cols], 1, which.max), levels = 1:length(topic_cols))
    
    # Handle Display_Name logic compatible with both Gene-as-word and SS-as-word
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
      geom_bar(stat="identity", position="stack", width = 0.85) + # Slightly thinner bars for better readability
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
        plot.margin = margin(20, 30, 20, 20) # Increase margin to prevent label clipping
      )
      
    # [Core Fix]: Increase width from 14 to 22 to eliminate horizontal compression
    ggsave(out_file, p, width=22, height=max(10, nrow(df_topics) * 0.5), limitsize = FALSE)
  }
  # ==========================================
  # 2. New: Top Words Distributions (SS-as-Word mapped to Trait)
  # ==========================================
  plot_mapped_top_words_gaussian <- function(gene_terms_mat, master_map, title_str, out_file, top_n = 5) {
    require(dplyr)
    require(ggplot2)
    require(stringr)

    # (A) Convert Gaussian weight matrix to long format and extract top_n
    df_long <- data.frame()
    for (i in 1:ncol(gene_terms_mat)) {
      topic_name <- colnames(gene_terms_mat)[i]
      if (!grepl("Topic", topic_name, ignore.case=TRUE)) {
        topic_name <- paste0("Topic ", i)
      }
      
      topic_weights <- gene_terms_mat[, i]
      # Sort descending and extract indices of top_n
      top_idx <- order(topic_weights, decreasing = TRUE)[1:top_n]
      
      temp_df <- data.frame(
        Topic = topic_name,
        Comorbidity = rownames(gene_terms_mat)[top_idx],
        Weight = topic_weights[top_idx],
        stringsAsFactors = FALSE
      )
      df_long <- rbind(df_long, temp_df)
    }

    # (B) Strip prefixes/suffixes and map to real names in master_map
    df_long$clean_id <- tolower(trimws(gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$|\\.txt$", "", df_long$Comorbidity)))
    
    if (!is.null(master_map)) {
      master_map$clean_filename <- tolower(trimws(gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$|\\.txt$", "", master_map$filename)))
      master_sub <- master_map[!duplicated(master_map$clean_filename), c("clean_filename", "trait")]
      
      df_mapped <- merge(df_long, master_sub, by.x = "clean_id", by.y = "clean_filename", all.x = TRUE)
      df_mapped$Display_Name <- ifelse(!is.na(df_mapped$trait) & trimws(df_mapped$trait) != "", 
                                       trimws(df_mapped$trait), df_mapped$Comorbidity)
    } else {
      df_mapped <- df_long
      df_mapped$Display_Name <- df_mapped$Comorbidity
    }
    
    # (C) Normalize weight processing
    plot_terms <- df_mapped %>%
      group_by(Topic) %>%
      arrange(desc(Weight)) %>% 
      slice_head(n = top_n) %>% 
      mutate(Normalized_Weight = Weight / max(Weight)) %>%
      ungroup() %>%
      mutate(unique_term_id = paste(Topic, Display_Name, sep = "___"))
    
    # (D) Smart layout logic (fixed 3 columns)
    n_topics <- length(unique(plot_terms$Topic))
    n_cols <- 3  
    n_rows <- ceiling(n_topics / n_cols)
    
    dynamic_width <- max(36, n_cols * 12) 
    dynamic_height <- max(12, n_rows * 6.0)
    
    # (E) Plot
    p_top_words <- ggplot(plot_terms, aes(x = reorder(unique_term_id, Normalized_Weight), y = Normalized_Weight, fill = Topic)) +
      geom_bar(stat = "identity") +
      facet_wrap(~Topic, scales = "free_y", ncol = n_cols) + 
      coord_flip() +
      theme_classic() +
      scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = c(0.01, 0)) +
      scale_x_discrete(labels = function(x) stringr::str_wrap(sub(".*___", "", x), width = 20)) +
      theme(
        text = element_text(size = 32),
        axis.text.y = element_text(size = 28, color = "black", face = "bold", lineheight = 0.9), 
        axis.text.x = element_text(size = 24, color = "black"),
        axis.title = element_text(size = 34, face = "bold", margin = margin(t = 15)),
        plot.title = element_text(size = 46, face = "bold", hjust = 0.5, margin = margin(b = 25)),
        strip.text = element_text(size = 32, face = "bold"),
        panel.spacing.x = unit(5.5, "lines"), 
        panel.spacing.y = unit(3.0, "lines"),
        plot.margin = margin(30, 30, 30, 30),
        legend.position = "none" # Hide legend to keep page clean
      ) +
      labs(title = title_str, x = "", y = "Relative Importance")
    
    ggsave(out_file, p_top_words, width = dynamic_width, height = dynamic_height, dpi = 300, limitsize = FALSE)
  }

  # Call plot function
  plot_topic_dist(gaussian_topic_matrix_std, "Gaussian Topic Distribution (Sorted)", paste0(prefix_gaussian, "_topic_dist_sorted.png"))
  
  # Call the Top 5 extraction function we just wrote
  plot_mapped_top_words_gaussian(out_gaussian$gene_terms, master_map, 
                                 "Top 5 Words per Topic (Gaussian Mapped & Normalized)", 
                                 paste0(prefix_gaussian, "_gaussian_top5_words_normalized.png"), 
                                 top_n = 5)

  cat("SS-as-Word Gaussian Analysis Complete.\n")
}