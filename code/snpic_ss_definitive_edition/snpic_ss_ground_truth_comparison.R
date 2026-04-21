############################################
# Compare observed similarity with ground truth (group-level)
############################################
# Build gene set-based similarity matrix from ground truth directory
# method = "overlap" uses |A∩B| / min(|A|,|B|); "jaccard" uses |A∩B| / |A∪B|
build_ground_truth_similarity_from_gt_dir <- function(ground_truth_dir,
                                                      snp_gene_map_file,
                                                      method = c("overlap","jaccard")) {
  method <- match.arg(method)
  
  # Read SNP → gene mapping table
  snp_gene_map <- data.table::fread(snp_gene_map_file, header = FALSE)
  colnames(snp_gene_map) <- c("snp", "gene")
  
  # Read all ground truth SNP files
  snp_files <- list.files(ground_truth_dir, pattern = "^snps_group.*\\.txt$", full.names = TRUE)
  if (length(snp_files) == 0) {
    stop("No snps_group*.txt found under: ", ground_truth_dir)
  }
  
  gene_sets <- list()
  for (f in snp_files) {
    # Row and column names need to be groupX_Y (remove snps_ prefix)
    group_name <- gsub("^snps_", "", tools::file_path_sans_ext(basename(f)))
    snps <- data.table::fread(f, header = FALSE)[[1]]
    genes <- unique(snp_gene_map$gene[match(snps, snp_gene_map$snp)])
    genes <- genes[!is.na(genes)]
    if (length(genes) == 0) {
      message("⚠️ No mapped genes for ", group_name, " (treated as empty set)")
      genes <- character(0)
    }
    gene_sets[[group_name]] <- genes
  }
  
  groups <- names(gene_sets)
  n <- length(groups)
  sim_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(groups, groups))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      A <- gene_sets[[i]]
      B <- gene_sets[[j]]
      inter <- length(intersect(A, B))
      if (method == "jaccard") {
        denom <- length(unique(c(A, B)))
      } else { # overlap / Simpson
        denom <- min(length(A), length(B))
      }
      sim_mat[i, j] <- ifelse(denom == 0, 0, inter / denom)
    }
  }
  return(sim_mat)
}

# Aggregate disease-level observed similarity to group-level (by label)
calculate_group_similarity <- function(observed_sim, disease_labels) {
  groups <- unique(disease_labels)
  group_sim <- matrix(0, nrow = length(groups), ncol = length(groups),
                      dimnames = list(groups, groups))
  for (g1 in groups) {
    for (g2 in groups) {
      members1 <- names(disease_labels)[disease_labels == g1]
      members2 <- names(disease_labels)[disease_labels == g2]
      vals <- observed_sim[members1, members2, drop = FALSE]
      group_sim[g1, g2] <- mean(vals)
    }
  }
  group_sim
}

# Calculate divergence between observed and ground truth (return difference matrix and mean absolute divergence)
calculate_group_divergence <- function(group_sim, ground_truth_sim) {
  common <- intersect(rownames(group_sim), rownames(ground_truth_sim))
  if (length(common) == 0) stop("No common group names between observed and ground truth matrices.")
  group_sim <- group_sim[common, common, drop = FALSE]
  ground_truth_sim <- ground_truth_sim[common, common, drop = FALSE]
  
  diff_mat <- group_sim - ground_truth_sim
  divergence <- mean(abs(diff_mat))
  list(divergence = divergence, diff_matrix = diff_mat)
}


############################################
# Calculate similarity for PCA and UMAP results
############################################
extract_coordinates <- function(umap_pca_results, method = c("pca", "umap")) {
  method <- match.arg(method)
  
  if (method == "pca") {
    coords_df <- umap_pca_results$pca_df
  } else {
    coords_df <- umap_pca_results$umap_df
  }
  
  # Extract numeric columns (PC1, PC2 or UMAP1, UMAP2)
  numeric_cols <- sapply(coords_df, is.numeric)
  coords <- coords_df[, numeric_cols, drop = FALSE]
  
  # Set row names as disease names
  if ("Disease" %in% colnames(coords_df)) {
    rownames(coords) <- coords_df$Disease
  } else if (is.null(rownames(coords))) {
    rownames(coords) <- paste0("sample_", 1:nrow(coords))
  }
  
  return(coords)
}

# Calculate similarity matrix for coordinate data
calculate_coordinate_similarity <- function(coords, method = "euclidean") {
  if (method == "euclidean") {
    dist_matrix <- dist(coords, method = "euclidean")
    # Convert distance to similarity (0-1 range)
    max_dist <- max(dist_matrix)
    sim <- 1 - (as.matrix(dist_matrix) / max_dist)
  } else if (method == "cosine") {
    sim <- proxy::simil(coords, method = "cosine")
    sim <- as.matrix(sim)
    # Map cosine similarity from [-1,1] to [0,1]
    sim <- (sim + 1) / 2
  } else {
    stop("Unsupported method. Choose from: euclidean, cosine")
  }
  
  return(sim)
}

# Analyze similarity for PCA and UMAP
analyze_dim_reduction_similarity <- function(umap_pca_results, disease_labels, ground_truth_sim) {
  results <- list()
  
  # Analyze PCA
  pca_coords <- extract_coordinates(umap_pca_results, "pca")
  pca_sim <- calculate_coordinate_similarity(pca_coords, "euclidean")
  pca_group_sim <- calculate_group_similarity(pca_sim, disease_labels)
  pca_div <- calculate_group_divergence(pca_group_sim, ground_truth_sim)
  results$pca <- list(similarity = pca_sim, group_similarity = pca_group_sim, divergence = pca_div)
  
  # Analyze UMAP
  umap_coords <- extract_coordinates(umap_pca_results, "umap")
  umap_sim <- calculate_coordinate_similarity(umap_coords, "euclidean")
  umap_group_sim <- calculate_group_similarity(umap_sim, disease_labels)
  umap_div <- calculate_group_divergence(umap_group_sim, ground_truth_sim)
  results$umap <- list(similarity = umap_sim, group_similarity = umap_group_sim, divergence = umap_div)
  
  return(results)
}
