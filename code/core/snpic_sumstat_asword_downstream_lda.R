# downstream_ss_lda.R
run_downstream_ss_lda <- function(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, seed, keep_all_traits) {
  cat("\n############################################\n")
  cat("Starting Downstream Analysis (SS-as-Word: Normal LDA)\n")
  cat("############################################\n")

  k_topic <- best_k
  prefix_downstream <- paste0(prefix, "_downstream_normal_k", k_topic)

  meta_df$Disease <- trimws(meta_df$Disease)
  passed_ids <- if (keep_all_traits) names(final_res$confidence_score) else names(final_res$confidence_score)[final_res$confidence_score >= best_thresh]
  passed_ids <- trimws(passed_ids)

  stable_list_df <- meta_df[meta_df$Disease %in% passed_ids, ]
  if (nrow(stable_list_df) < 2) stop("Error: Not enough stable diseases (<2) for downstream.")

  target_diseases <- stable_list_df$Disease
  rownames(input_mat) <- trimws(rownames(input_mat))
  common_diseases <- intersect(rownames(input_mat), target_diseases)

  # For SS-as-word, the matrix is symmetric; trim rows and columns simultaneously
  matrix_result <- input_mat[common_diseases, common_diseases, drop = FALSE]
  valid_idx <- apply(matrix_result, 2, function(x) length(unique(x)) > 1)
  matrix_result <- matrix_result[valid_idx, valid_idx, drop = FALSE]
  cat(sprintf(">> Filtered SS Matrix Dimensions: %d Diseases x %d Diseases (Words)\n", nrow(matrix_result), ncol(matrix_result)))

  # Run SNPic (LDA)
  # set.seed(seed)
  # results_umap_pca <- plot_umap_pca_from_matrix(matrix_result, prefix_downstream, master_map)
  set.seed(seed)
  results_snpic <- run_lda_analysis(matrix_result, k = k_topic, prefix = prefix_downstream, master_map = master_map)
  
  snpic_topic_matrix <- results_snpic$topic_matrix
  snpic_top_word <- results_snpic$top_word

  # Similarity Analysis
  sim_res_lda <- analyze_disease_similarity(snpic_topic_matrix, model_type = "LDA", mode = "separate", edge_threshold = 0.70, seed = seed)
  if(!is.null(sim_res_lda$correlation_plots$network)) ggsave(paste0(prefix_downstream, "_similarity_network_lda.png"), sim_res_lda$correlation_plots$network, width=14, height=10)

  # Topic Distributions (SS-as-word custom version: larger font, narrower layout)
  plot_topic_dist <- function(topic_mat, title_str, out_file) {
    non_topic_cols <- c("Disease", "disease_clean", "label_map", "label", "label2", "Disease_file", "Disease_Unique", "MaxTopic")
    topic_cols <- setdiff(colnames(topic_mat), non_topic_cols)
    
    cols_to_keep <- c(topic_cols, "Disease", if ("label2" %in% colnames(topic_mat)) "label2" else NULL)
    df_topics <- topic_mat[, cols_to_keep]
    df_topics$MaxTopic <- factor(apply(df_topics[, topic_cols], 1, which.max), levels = 1:length(topic_cols))
    
    # This handles both Gene-as-word and SS-as-word Display_Name logic
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
      geom_bar(stat="identity", position="stack", width = 0.85) + # Slightly thinner bars for better breathing room
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
        plot.margin = margin(20, 30, 20, 20) # Add margin around plot to prevent label clipping
      )
      
    # [Core Fix]: Increased width from 14 to 22 to completely resolve horizontal compression
    ggsave(out_file, p, width=22, height=max(10, nrow(df_topics) * 0.5), limitsize = FALSE)
  }  
  # Top Words Mapping (SS-as-word specific)
  plot_mapped_top_words <- function(top_word_df, master_map, title_str, out_file) {
    top_word_df$clean_id <- tolower(trimws(gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$|\\.txt$", "", top_word_df$Comorbidity)))
    if (!is.null(master_map)) {
      master_map$clean_filename <- tolower(trimws(gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$|\\.txt$", "", master_map$filename)))
      master_sub <- master_map[!duplicated(master_map$clean_filename), c("clean_filename", "trait")]
      df_mapped <- merge(top_word_df, master_sub, by.x = "clean_id", by.y = "clean_filename", all.x = TRUE)
      df_mapped$Display_Name <- ifelse(!is.na(df_mapped$trait) & trimws(df_mapped$trait) != "", trimws(df_mapped$trait), df_mapped$Comorbidity)
    } else {
      df_mapped <- top_word_df
      df_mapped$Display_Name <- df_mapped$Comorbidity
    }
    
    plot_terms <- df_mapped %>% group_by(Topic) %>% arrange(desc(Weight)) %>% slice_head(n = 5) %>% 
      mutate(Normalized_Weight = Weight / max(Weight)) %>% ungroup() %>% mutate(unique_term_id = paste(Topic, Display_Name, sep = "___"))
    
    n_topics <- length(unique(plot_terms$Topic)); n_cols <- 3; n_rows <- ceiling(n_topics / n_cols)
    p_top_words <- ggplot(plot_terms, aes(x = reorder(unique_term_id, Normalized_Weight), y = Normalized_Weight, fill = Topic)) +
      geom_bar(stat = "identity") + facet_wrap(~Topic, scales = "free_y", ncol = n_cols) + coord_flip() + theme_classic() +
      scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1), expand = c(0.01, 0)) +
      scale_x_discrete(labels = function(x) stringr::str_wrap(sub(".*___", "", x), width = 20)) +
      theme(
        text = element_text(size = 32), axis.text.y = element_text(size = 28, color = "black", face = "bold", lineheight = 0.9), 
        axis.text.x = element_text(size = 24, color = "black"), axis.title = element_text(size = 34, face = "bold", margin = margin(t = 15)),
        plot.title = element_text(size = 46, face = "bold", hjust = 0.5, margin = margin(b = 25)),
        strip.text = element_text(size = 32, face = "bold"), panel.spacing.x = unit(5.5, "lines"), panel.spacing.y = unit(3.0, "lines"), plot.margin = margin(30, 30, 30, 30)
      ) + labs(title = title_str, x = "", y = "Relative Importance")
    ggsave(out_file, p_top_words, width = max(40, n_cols * 12), height = max(12, n_rows * 6.0), dpi = 300, limitsize = FALSE)
  }

  plot_topic_dist(snpic_topic_matrix, "SNPic LDA Topic Distribution (Sorted)", paste0(prefix_downstream, "_topic_dist_sorted.png"))
  plot_mapped_top_words(snpic_top_word, master_map, "Top 5 Words per Topic (Mapped & Normalized)", paste0(prefix_downstream, "_snpic_top5_words_normalized.png"))

  cat("SS-as-Word Normal LDA Analysis Complete.\n")
}