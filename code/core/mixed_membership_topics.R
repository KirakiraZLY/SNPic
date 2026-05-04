library(maptpx)
library(reshape2)
library(ggplot2)
library(dplyr)
library(umap)
set.seed(123)
run_mixed_membership <- function(matrix_result, k = 5, custom_order = NULL, master_map = NULL, standardize_by_row = TRUE,
                                 # Tunable hyperparameters
                                 tol = 0.01,         # Convergence threshold
                                 shape = 0.01,       # Prior shape parameter
                                 tmax = 2000,        # Maximum iterations
                                 verb = 1,           # Verbose level
                                 # New adjustable parameters
                                 init_type = "random", # Initialization type
                                 prior_weight = 0.1) { # Prior weight
  
  # Data standardization
  if (standardize_by_row) {
    # Standardize by row
    dtm_scaled <- as.data.frame(t(scale(t(matrix_result))))
  } else {
    # Standardize by column
    dtm_scaled <- as.data.frame(scale(matrix_result))
  }

  # Ensure non-negative
  dtm_scaled <- dtm_scaled - min(dtm_scaled, na.rm = TRUE)

  # Convert to matrix, ensure non-negative
  dtm <- as.matrix(dtm_scaled)
  dtm[dtm < 0] <- 0
  # Remove all-zero rows
  dtm <- dtm[rowSums(dtm != 0) > 0, ]
  
  # Convert to DocumentTermMatrix
  dtm_dtm <- as.DocumentTermMatrix(
    slam::as.simple_triplet_matrix(dtm),
    weighting = weightTf
  )
  
  # Run topics()
  mmt_model <- maptpx::topics(
    dtm_dtm,
    K = k,
    tol = tol,
    verb = verb,
    tmax = tmax,
    shape = shape,
    # Possible additional parameters
    initopics = NULL,  # Can specify initial topics
    bf = FALSE         # Whether to use BFGS optimization
  )
  
  # Step 3: Extract posterior topic distribution per disease (theta)
  topics_matrix <- mmt_model$omega  # rows: diseases, cols: topics
  rownames(topics_matrix) <- rownames(dtm)
  colnames(topics_matrix) <- paste0("Topic", 1:k)
  
  # Step 4: Melt for stacked barplot
  df_topics <- as.data.frame(topics_matrix)
  df_topics$Disease <- rownames(df_topics)
  df_melted <- reshape2::melt(df_topics, id.vars = "Disease", variable.name = "Topic", value.name = "Probability")
  
  # Modify first plot - rotate 90 degrees
  gg1 <- ggplot(df_melted, aes(x = Probability, y = Disease, fill = Topic)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_minimal() +
    labs(title = "Disease Topic distribution across Mixed-Membership Model",
         x = "Topic Probability", y = "Disease") +
    theme(axis.text.y = element_text(size = 14),
          legend.position = "right")
  # print(gg1)
  
  # Step 5: MaxTopic processing and second plot
  df_topics$MaxTopic <- apply(topics_matrix, 1, which.max)
  
  # Use custom order (default sort if not provided)
  if (is.null(custom_order)) {
    custom_order <- sort(unique(df_topics$MaxTopic))
  }
  
  df_topics$MaxTopic <- factor(df_topics$MaxTopic,
                               levels = custom_order,
                               labels = paste0("Topic", custom_order))
  
  # Melt only topic probability columns
  topic_cols <- grep("^Topic\\d+$", names(df_topics), value = TRUE) # Get all topic columns
  
  df_melted <- df_topics %>%
    pivot_longer(
      cols = all_of(topic_cols),
      names_to = "Topic",
      values_to = "Probability"
    )
  
  # Modify second plot - rotate 90 degrees
  gg2 <- ggplot(df_melted, aes(x = Probability, y = Disease, fill = Topic)) +
    geom_col(position = "stack") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Disease Topic distribution across Mixed-Membership Model",
      subtitle = "Sorted by max topic probability",
      x = "Topic Probability",
      y = "Disease"
    ) +
    theme(
      axis.text.y = element_text(size = 14),
      panel.grid.major.y = element_blank(),  # Changed to major.y
      plot.margin = margin(10, 10, 10, 10),
      legend.position = "right"  # Legend on the right
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  # Add facets when too many diseases - adjust facet direction
  if (n_distinct(df_melted$Disease) > 30) {
    gg2 <- gg2 + 
      facet_grid(MaxTopic ~ ., scales = "free_y", space = "free_y") +  # Changed to vertical facets
      theme(
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 8)
      )
  }
  
  # print(gg2)
  topic_matrix2 <- as.data.frame(topics_matrix)
  topic_matrix2$Disease <- rownames(topics_matrix)
  topic_matrix2 <- resolve_labels(topic_matrix2, master_map)
  
  return(list(
    model = mmt_model,
    topic_matrix = topics_matrix,
    topic_matrix2 = topic_matrix2,
    gene_terms = mmt_model$theta,  # terms per topic
    topic_df = df_melted,
    standardized_matrix = dtm_scaled  # Return standardized matrix
  ))
}

# Unified and clean resolve_labels function
resolve_labels <- function(df, master_map = NULL) {
  if (is.null(master_map)) {
    df$label <- "Other"
    df$label2 <- "Unknown"
    return(df)
  }
  
  # Clean disease names for matching
  if (!"disease_clean" %in% colnames(df)) {
    df$disease_clean <- gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$", "", df$Disease)
  }
  
  # Prevent data inflation from duplicate matching
  master_sub <- master_map[!duplicated(master_map$clean_filename), c("clean_filename", "trait", "label1", "label2")]
  
  # Merge
  df_merged <- merge(df, master_sub, by.x = "disease_clean", by.y = "clean_filename", all.x = TRUE)
  
  # Override: use trait if available, label1 if available, otherwise defaults
  df_merged$Disease <- ifelse(!is.na(df_merged$trait) & df_merged$trait != "", df_merged$trait, df_merged$Disease)
  df_merged$label <- ifelse(!is.na(df_merged$label1) & df_merged$label1 != "", df_merged$label1, "Other")
  df_merged$label2 <- ifelse(!is.na(df_merged$label2) & df_merged$label2 != "", df_merged$label2, "Unknown")
  
  # Remove auxiliary columns
  df_merged$disease_clean <- NULL
  df_merged$trait <- NULL
  df_merged$label1 <- NULL
  
  return(df_merged)
}

plot_top_genes <- function(gene_terms, top_n = 10) {
  gene_df <- as.data.frame(gene_terms)
  colnames(gene_df) <- paste0("Topic", 1:ncol(gene_df))
  gene_df$Gene <- rownames(gene_df)
  
  gene_long <- reshape2::melt(gene_df, id.vars = "Gene", variable.name = "Topic", value.name = "Weight")
  
  # Select top_n genes per topic
  top_genes <- gene_long %>%
    group_by(Topic) %>%
    top_n(top_n, Weight) %>%
    ungroup()
  
  gg <- ggplot(top_genes, aes(x = reorder(Gene, Weight), y = Weight, fill = Topic)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    facet_wrap(~ Topic, scales = "free") +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme_classic() +
    theme(
      text = element_text(size = 12),           # All text font size is 12
      axis.text = element_text(size = 12),      # Axis text font size is 12
      axis.title = element_text(size = 12),     # Axis title font size is 12
      plot.title = element_text(size = 12),     # Plot title font size is 12
      strip.text = element_text(size = 12)      # Facet strip text font size is 12
    ) +
    labs(
      title = paste("Top", top_n, "genes for each topic"),
      x = "Gene",
      y = "Weight in Topic"
    )
  print(gg)
  return(gg)
}

plot_umap_topics <- function(topic_matrix, 
                             prefix = NULL,  
                             map_detail_func = NULL) {
  # topic_matrix <- out$topic_matrix
  ## UMAP
  # Run UMAP on topic_matrix
  umap_config <- umap::umap.defaults
  umap_config$n_neighbors <- 5
  umap_config$random_state <- 123
  umap_res <- umap(topic_matrix, config = umap_config)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Disease <- rownames(topic_matrix)
  
  # label
  umap_df <- resolve_labels(umap_df, master_map)
  
  # Get color configuration
  custom_colors <- get_network_colors()
  unique_labels <- unique(umap_df$label)
  n_labels <- length(unique_labels)
  
  # visualization
  p_snpic_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = label)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Disease), size = 3, max.overlaps = Inf) +
    theme_classic() +
    ggtitle("SS as word, SNPic-SS + UMAP Clustering, by label") +
    scale_color_manual(
      values = if(n_labels <= length(custom_colors)) {
        setNames(custom_colors[1:n_labels], unique_labels)
      } else {
        NULL
      }
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )
  # print(p_snpic_umap)
  # # optional: save plot
  # if (!is.null(prefix)) {
  #   ggsave(paste0(prefix, "_SSword_snpic_umap.png"), p_snpic_umap, width = 10, height = 8)
  # }
  return(umap_df)
}

plot_pca_topics <- function(topic_matrix, prefix = NULL, master_map = NULL) {  
  pca <- prcomp(topic_matrix, center = TRUE, scale. = TRUE)
  df_pca <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Disease = rownames(topic_matrix),
    MaxTopic = apply(topic_matrix, 1, which.max)
  )
  
  # Directly use resolve_labels to set labels
  df_pca <- resolve_labels(df_pca, master_map)
  
  # Get color configuration
  custom_colors <- get_network_colors()
  unique_labels <- unique(df_pca$label)
  n_labels <- length(unique_labels)
  
  p_snpic_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = label)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = Disease), size = 3, max.overlaps = 20) +
    theme_minimal() +
    labs(
      title = "SNPic + PCA projection of diseases based on topic distribution",
      x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)")
    ) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    scale_color_manual(
      values = if(n_labels <= length(custom_colors)) {
        setNames(custom_colors[1:n_labels], unique_labels)
      } else {
        NULL
      }
    )
  # print(p_snpic_pca)
  # if (!is.null(prefix)) {
  #   ggsave(paste0(prefix, "_SSword_snpic_pca.png"), p_snpic_pca, width = 10, height = 8)
  # }
  return(list(pca_df = df_pca, pca_obj = pca))
}