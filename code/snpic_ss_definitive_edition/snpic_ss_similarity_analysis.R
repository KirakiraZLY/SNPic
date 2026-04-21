############################################
# Similarity calculation and visualization functions (Modified for Duplicates & Legend)
############################################

# 2. 计算相似度 (不变)
calculate_disease_similarity <- function(topic_numeric, method = "correlation") {
  if(!is.matrix(topic_numeric)) {
    topic_numeric <- as.matrix(topic_numeric)
  }
  topic_numeric[is.na(topic_numeric)] <- 0
  row_sums <- rowSums(topic_numeric)
  if(any(row_sums == 0)) {
    warning(paste("Removing", sum(row_sums == 0), "rows with all zero values."))
    topic_numeric <- topic_numeric[row_sums > 0, , drop = FALSE]
  }
  
  if(method == "hellinger") {
    topic_prob <- sweep(topic_numeric, 1, rowSums(topic_numeric), FUN = "/")
    dist_matrix <- dist(sqrt(topic_prob), method = "euclidean")
    sim <- 1 - (as.matrix(dist_matrix) / sqrt(2))
  } else if(method == "js") {
    P <- sweep(topic_numeric, 1, rowSums(topic_numeric), FUN = "/")
    n <- nrow(P)
    sim <- matrix(0, nrow = n, ncol = n)
    rownames(sim) <- rownames(topic_numeric)
    colnames(sim) <- rownames(topic_numeric)
    js_dist <- function(p, q) {
      m <- 0.5 * (p + q)
      kl_pm <- sum(p * log((p + 1e-12) / (m + 1e-12)))
      kl_qm <- sum(q * log((q + 1e-12) / (m + 1e-12)))
      return(0.5 * (kl_pm + kl_qm))
    }
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        val <- 1 - sqrt(js_dist(P[i, ], P[j, ]))
        sim[i, j] <- val
        sim[j, i] <- val
      }
    }
    diag(sim) <- 1
  } else if(method == "correlation") {
    sd_vec <- apply(topic_numeric, 1, sd)
    if(any(sd_vec == 0)) topic_numeric <- topic_numeric[sd_vec > 0, ]
    corr <- cor(t(topic_numeric))
    corr[is.na(corr)] <- 0
    sim <- (corr + 1) / 2
  } else {
    stop("Unsupported method.")
  }
  return(sim)
}

# 3. 热图可视化 (不变)
visualize_topic_heatmap <- function(topic_numeric, disease_labels = NULL, model_type = "Unknown", method = "correlation") {
  require(ggplot2)
  require(reshape2)
  heatmap_data <- as.data.frame(topic_numeric)
  heatmap_data$UniqueID <- rownames(heatmap_data)
  
  if(!is.null(disease_labels)) {
    heatmap_data$Group <- disease_labels[match(heatmap_data$UniqueID, names(disease_labels))]
    heatmap_data$Group[is.na(heatmap_data$Group)] <- "Unknown"
  } else {
    heatmap_data$Group <- sapply(heatmap_data$UniqueID, function(x) {
      clean_name <- gsub("\\.\\d+$", "", x)
      parts <- strsplit(clean_name, "_")[[1]]
      if(length(parts) >= 2) parts[1] else "Unknown"
    })
  }
  
  heatmap_data <- heatmap_data[order(heatmap_data$Group, heatmap_data$UniqueID), ]
  heatmap_data$UniqueID <- factor(heatmap_data$UniqueID, levels = unique(heatmap_data$UniqueID))
  
  melted_data <- reshape2::melt(heatmap_data, id.vars = c("UniqueID", "Group"))
  
  ggplot(melted_data, aes(x = variable, y = UniqueID, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5, name = "Prob/Value") +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 6),
          strip.text.y = element_text(angle = 0, size = 8, face = "bold"),
          panel.spacing = unit(0.1, "lines")) +
    labs(title = paste("Topic Distribution Heatmap -", model_type, "Model"),
         subtitle = paste("Similarity method:", method),
         x = "Topic", y = "Trait (Unique ID)")
}

# 4. 网络图 (Shape Only) - 添加 Edge Weight 图例
visualize_network_shape_only <- function(similarity_matrix, disease_labels = NULL, model_type = "Unknown", method = "correlation", topic_matrix = NULL, edge_threshold_percentile = 0.70, seed = 123) {
  require(igraph)
  require(ggrepel)
  require(ggplot2)
  
  threshold_val <- quantile(similarity_matrix[upper.tri(similarity_matrix)], edge_threshold_percentile)
  adj <- similarity_matrix
  adj[adj < threshold_val] <- 0
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  if(vcount(g) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data available"))
  }
  
  node_ids <- V(g)$name
  groups <- rep("Unknown", length(node_ids))
  display_names <- node_ids
  
  if(!is.null(topic_matrix) && "UniqueID" %in% colnames(topic_matrix)) {
    idx <- match(node_ids, topic_matrix$UniqueID)
    
    if("Trait" %in% colnames(topic_matrix)) {
      display_names <- topic_matrix$Disease[idx]
    }
    if("label2" %in% colnames(topic_matrix)) {
      matched_labels <- topic_matrix$label2[idx]
      groups[!is.na(matched_labels)] <- matched_labels[!is.na(matched_labels)]
    }
  }
  
  V(g)$group <- groups
  
  set.seed(seed)
  layout <- layout_with_fr(g)
  
  node_pos <- data.frame(
    x = layout[, 1],
    y = layout[, 2],
    name = V(g)$name,
    label = display_names,
    group = V(g)$group
  )
  
  edges <- as.data.frame(get.edgelist(g))
  has_edges <- nrow(edges) > 0 
  
  if(has_edges) {
    colnames(edges) <- c("from", "to")
    edges$weight <- E(g)$weight
    edges <- merge(edges, node_pos, by.x = "from", by.y = "name")
    edges <- merge(edges, node_pos, by.x = "to", by.y = "name", suffixes = c(".from", ".to"))
  }
  
  unique_groups <- unique(node_pos$group)
  custom_colors <- get_network_colors()
  n_groups <- length(unique_groups)
  my_shapes <- rep(c(16, 17, 15, 18, 4, 7, 3, 10, 11, 12, 13, 14), length.out = n_groups)
  my_colors <- rep(custom_colors, length.out = n_groups)
  
  shape_values <- setNames(my_shapes, unique_groups)
  color_values <- setNames(my_colors, unique_groups)
  
  p <- ggplot() +
    {if(has_edges) 
      geom_segment(data = edges,
                   aes(x = x.from, y = y.from, xend = x.to, yend = y.to, 
                       alpha = weight, linewidth = weight), 
                   color = "grey60") 
    } +
    # 【修改】：点的大小从 4 改大到 6
    geom_point(data = node_pos,
               aes(x = x, y = y, color = group, shape = group),
               size = 6, alpha = 0.9) +
    geom_text_repel(data = node_pos, 
                    aes(x = x, y = y, label = label), 
                    size = 3.5,
                    max.overlaps = 50, 
                    bg.color = "white", 
                    bg.r = 0.15) +
    scale_alpha(range = c(0.2, 0.8), guide = "none") + 
    scale_linewidth_continuous(name = "Edge Weight", range = c(0.2, 2.0)) +  
    scale_color_manual(name = "Data Source", values = color_values) +
    scale_shape_manual(name = "Data Source", values = shape_values) +
    theme_void() +
    theme(
      text = element_text(size = 14),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(title = paste("Trait Similarity Network -", model_type),
         subtitle = paste("Method:", method, "| Edge threshold percentile:", edge_threshold_percentile * 100, "%"))
  
  return(p)
}

# 5. 网络图 (Separate) - 添加 Edge Weight 图例，点放大
visualize_network_separate <- function(similarity_matrix, disease_labels = NULL, model_type = "Unknown", method = "correlation", topic_matrix = NULL, edge_threshold_percentile = 0.70, seed = 123) {
  require(igraph)
  require(ggrepel)
  require(ggplot2)
  
  threshold_val <- quantile(similarity_matrix[upper.tri(similarity_matrix)], edge_threshold_percentile)
  adj <- similarity_matrix
  adj[adj < threshold_val] <- 0
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  if(vcount(g) == 0) return(ggplot() + theme_void() + labs(title = "No data available"))
  
  node_ids <- V(g)$name
  color_groups <- rep("Unknown", length(node_ids))
  shape_groups <- rep("Unknown", length(node_ids))
  display_names <- node_ids
  
  if(!is.null(topic_matrix) && "UniqueID" %in% colnames(topic_matrix)) {
    idx <- match(node_ids, topic_matrix$UniqueID)
    
    if("Trait" %in% colnames(topic_matrix)) {
      display_names <- topic_matrix$Disease[idx] 
    }
    
    if("label" %in% colnames(topic_matrix)) {
      vals <- topic_matrix$label[idx]
      color_groups[!is.na(vals)] <- vals[!is.na(vals)]
    }
    if("label2" %in% colnames(topic_matrix)) {
      vals <- topic_matrix$label2[idx]
      shape_groups[!is.na(vals)] <- vals[!is.na(vals)]
    }
  }
  
  V(g)$color_group <- color_groups
  V(g)$shape_group <- shape_groups
  
  set.seed(seed) 
  layout <- layout_with_fr(g)
  
  node_pos <- data.frame(
    x = layout[, 1],
    y = layout[, 2],
    name = V(g)$name,
    label = display_names,
    color_group = V(g)$color_group,
    shape_group = V(g)$shape_group
  )
  
  edges <- as.data.frame(get.edgelist(g))
  has_edges <- nrow(edges) > 0
  
  if(has_edges) {
    colnames(edges) <- c("from", "to")
    edges$weight <- E(g)$weight
    edges <- merge(edges, node_pos, by.x = "from", by.y = "name")
    edges <- merge(edges, node_pos, by.x = "to", by.y = "name", suffixes = c(".from", ".to"))
  }
  
  unique_color_groups <- unique(node_pos$color_group)
  unique_shape_groups <- unique(node_pos$shape_group)
  custom_colors <- get_network_colors()
  
  color_values <- setNames(rep(custom_colors, length.out = length(unique_color_groups)), unique_color_groups)
  shape_values <- setNames(rep(c(16, 17, 15, 18, 4, 7, 3, 10, 11, 12, 13, 14), length.out = length(unique_shape_groups)), unique_shape_groups)
  
  p <- ggplot() +
    {if(has_edges) 
      geom_segment(data = edges,
                   aes(x = x.from, y = y.from, xend = x.to, yend = y.to, 
                       alpha = weight, linewidth = weight),
                   color = "grey60")
    } +
    # 【修改】：点的大小从 4 改大到 6
    geom_point(data = node_pos,
               aes(x = x, y = y, color = color_group, shape = shape_group),
               size = 6, alpha = 0.9) +
    geom_text_repel(data = node_pos, 
                    aes(x = x, y = y, label = label), 
                    size = 3.5,
                    max.overlaps = 50,
                    bg.color = "white",
                    bg.r = 0.15) +
    scale_alpha(range = c(0.2, 0.8), guide = "none") +
    scale_linewidth_continuous(name = "Edge Weight", range = c(0.2, 2.0)) +
    scale_color_manual(name = "Trait Category", values = color_values) +
    scale_shape_manual(name = "Data Source", values = shape_values) +
    theme_void() +
    theme(
      text = element_text(size = 14),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = paste("Trait Similarity Network -", model_type),
         subtitle = paste("Method:", method, "| Edge threshold percentile:", edge_threshold_percentile * 100, "%"))
  
  return(p)
}

# 6. 主分析函数 (不变)
analyze_disease_similarity <- function(topic_matrix, methods = c("hellinger", "correlation"), 
                                       model_type = "Unknown", mode = "shape_only", 
                                       edge_threshold = 0.70, seed = 123,
                                       text_size = 6) { 
  
  require(ggplot2)
  results <- list()
  
  topic_matrix$UniqueID <- make.unique(as.character(topic_matrix$Disease))
  
  numeric_cols <- sapply(topic_matrix, is.numeric)
  topic_numeric <- topic_matrix[, numeric_cols, drop = FALSE]
  
  rownames(topic_numeric) <- topic_matrix$UniqueID
  
  topic_numeric[is.na(topic_numeric)] <- 0
  topic_numeric[topic_numeric < 0] <- 0
  
  if("label" %in% colnames(topic_matrix)) {
    disease_labels <- topic_matrix$label
    names(disease_labels) <- topic_matrix$UniqueID 
  } else {
    disease_labels <- NULL
  }
  
  for(method in methods) {
    message(paste("Calculating similarity using", method, "for", model_type, "model..."))
    sim <- calculate_disease_similarity(topic_numeric, method)
    results[[method]] <- sim
    
    message(paste("Visualizing results for", method, "..."))
    heatmap <- visualize_topic_heatmap(topic_numeric, disease_labels, model_type, method)
    
    if(mode == "shape_only") {
      network <- visualize_network_shape_only(sim, disease_labels, model_type, method, topic_matrix, 
                                              edge_threshold_percentile = edge_threshold, seed = seed)
    } else if(mode == "separate") {
      network <- visualize_network_separate(sim, disease_labels, model_type, method, topic_matrix, 
                                            edge_threshold_percentile = edge_threshold, seed = seed)
    } else {
      stop("Unsupported mode. Choose from: shape_only, separate")
    }
    
    if (!is.null(network)) {
      for(i in seq_along(network$layers)) {
        if(inherits(network$layers[[i]]$geom, "GeomTextRepel") || 
           inherits(network$layers[[i]]$geom, "GeomText")) {
          network$layers[[i]]$aes_params$size <- text_size
        }
      }
    }
    
    print(heatmap)
    print(network)
    results[[paste0(method, "_plots")]] <- list(heatmap = heatmap, network = network)
  }
  return(results)
}