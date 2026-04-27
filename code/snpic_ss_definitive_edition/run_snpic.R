#!/usr/bin/env Rscript

############################################
# 0. Parse Command Line Arguments
############################################
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
}
library(optparse)

option_list = list(
  make_option(c("-i", "--input_folder"), type="character", default=NULL, 
              help="[REQUIRED] Input folder containing disease/SNP list files", metavar="character"),
  make_option(c("-m", "--master_map"), type="character", default=NULL, 
              help="[REQUIRED] Path to master_disease_mapping.csv, including [filename,trait,label1,label2]", metavar="character"),
  make_option(c("-o", "--out_prefix"), type="character", default=NULL, 
              help="[REQUIRED] Output prefix (including path) for results", metavar="character"),
  make_option(c("-s", "--snp_gene_map"), type="character", default=NULL, 
              help="[OPTIONAL] Path to the SNP-gene projection map file (REQUIRED ONLY IF --mode is 'gene')", metavar="character"),
  
  make_option(c("--k_min"), type="integer", default=5, 
              help="Minimum number of topics in K range [default= %default]", metavar="integer"),
  make_option(c("--k_max"), type="integer", default=10, 
              help="Maximum number of topics in K range [default= %default]", metavar="integer"),
  make_option(c("-b", "--n_bootstrap"), type="integer", default=50, 
              help="Number of bootstrap iterations per K [default= %default]", metavar="integer"),
  make_option(c("-c", "--cores"), type="integer", default=NULL, 
              help="Number of CPU cores to use. Default is auto-detect.", metavar="integer"),
  make_option(c("--seed"), type="integer", default=123, 
              help="Random seed for reproducibility [default= %default]", metavar="integer"),
  
  make_option(c("--keep_all_traits"), action="store_true", default=FALSE, 
              help="If set, skip the stability filter and keep all traits for downstream analysis. [default= FALSE]"),
  
  make_option(c("--k_only"), type="integer", default=NULL, 
              help="[OPTIONAL] Fast mode. Skip bootstrap and run exactly this many topics. Automatically implies --keep_all_traits. [default= %default]", metavar="integer"),
  
  make_option(c("--mode"), type="character", default="gene", 
              help="Base matrix construction mode: 'gene' (Gene-as-word) or 'ss' (Sumstat-as-word). [default= %default]"),
  make_option(c("--model"), type="character", default="lda", 
              help="Downstream model to execute: 'lda', 'gaussian', or 'both'. [default= %default]")
)

opt_parser = OptionParser(option_list=option_list, description="Global SNPic Pipeline (Gene/SS x LDA/Gaussian Dispatcher)")
opt = parse_args(opt_parser)

# Base validations (Removed "dependency")
required_args <- c("input_folder", "master_map", "out_prefix")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(sprintf("Argument --%s is missing and is required.", arg), call.=FALSE)
  }
}

if (!(opt$mode %in% c("gene", "ss"))) stop("Invalid --mode argument. Choose 'gene' or 'ss'.")
if (!(opt$model %in% c("lda", "gaussian", "both"))) stop("Invalid --model argument. Choose 'lda', 'gaussian', or 'both'.")

if (opt$mode == "gene" && is.null(opt$snp_gene_map)) {
  stop("Error: --snp_gene_map is REQUIRED when --mode is 'gene'.")
}

############################################
# 1. Load dependencies & Setup (Auto-detect Path)
############################################
# 自动获取当前主脚本所在目录
args_info <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", args_info, value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(sub("--file=", "", file_arg[1]))
} else {
  script_dir <- getwd() # 如果不是通过 Rscript 运行，默认使用当前工作目录
}

dependency_file <- file.path(script_dir, "snpic_ss_dependency.R")
cat(sprintf("Loading dependencies dynamically from: %s\n", dependency_file))

if (file.exists(dependency_file)) {
  source(dependency_file)
} else {
  stop(sprintf("Error: Dependency file 'snpic_ss_dependency.R' not found in %s. Please ensure it is in the same directory as run_snpic.R.", script_dir))
}
set.seed(opt$seed)

############################################
# 2-4. Init Setup & Base Matrix Construction
############################################
folder1 <- opt$input_folder
prefix <- opt$out_prefix
output_dir <- dirname(prefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

k_range <- c(opt$k_min:opt$k_max)
n_bootstrap <- opt$n_bootstrap
master_map_file <- opt$master_map

# === 【核心路由一】：矩阵生成模式 ===
if (opt$mode == "gene") {
  cat("Constructing Gene-as-Word Matrix...\n")
  matrix_result <- run_gene_as_word_analysis(folder = folder1, snp_gene_map_file = opt$snp_gene_map, prefix = NULL, top_n_genes = 99999)
} else {
  cat("\nConstructing Shared SNP Matrix (SS-as-Word)...\n")
  matrix_result <- calculate_shared_snp_matrix(folder1 = folder1, folder2 = folder1, plot_heatmap = FALSE, plot_dendrogram = FALSE, prefix = NULL)
}

# 清洗与格式化
matrix_result <- matrix_result[apply(matrix_result, 1, function(x) !all(x == 0)), ]
matrix_result <- matrix_result[, apply(matrix_result, 2, function(x) !all(x == 0))]
input_mat <- as.matrix(matrix_result)

cat("Loading Master Mapping Table...\n")
if (file.exists(master_map_file)) {
  master_map <- read.csv(master_map_file, stringsAsFactors = FALSE)
  colnames(master_map)[colnames(master_map) == "disease"] <- "trait"
  master_map$clean_filename <- gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$", "", trimws(master_map$filename))
} else {
  master_map <- NULL
}

if (!exists("get_network_colors")) {
  get_network_colors <- function() c(
    "#C65D57", "#3A5A82", "#D4A853", "#D0E3A4", "#5D4037", "#664673", 
    "#9C3E3A", "#587CA8", "#C69B4B", "#85A143", "#795548", "#8A6796",
    "#D97A76", "#7395BF", "#DFB769", "#9FB860", "#8D6E63", "#9E7FA8",
    "#E59A96", "#91B0D1", "#EAC782", "#607A2A", "#A1887F", "#BCA3C4",
    "#F0B8B5", "#B6CEE6", "#F3D69C", "#B8CD81", "#BCAAA4", "#D4C2D9",
    "#A84F4B", "#476994", "#A8833C", "#708D33", "#4E342E", "#765582"
  )
}

run_lda_bootstrap_parallel <- function(dtm_matrix, k, n_iter = 50) {
  diseases <- rownames(dtm_matrix)
  n_disease <- length(diseases)
  n_cores <- if (!is.null(opt$cores)) opt$cores else max(1, min(parallel::detectCores() - 2, 20))
  cat(sprintf("  -> Running Bootstrap (k=%d, n=%d) on %d cores...\n", k, n_iter, n_cores))
  
  run_single <- function(i) {
    set.seed(opt$seed + i * 1000) 
    resampled_matrix <- dtm_matrix * 0 
    for (r in 1:nrow(dtm_matrix)) {
      counts <- dtm_matrix[r, ]
      total <- sum(counts)
      if (total > 0) resampled_matrix[r, ] <- as.vector(rmultinom(1, total, counts/total))
    }
    valid <- rowSums(resampled_matrix) > 0
    if(sum(valid) < 2) return(NULL)
    
    curr_dtm <- tm::as.DocumentTermMatrix(resampled_matrix[valid, ], weighting = tm::weightTf)
    tryCatch({
      mod <- topicmodels::LDA(curr_dtm, k = k, method = "Gibbs", control = list(seed = 1234+i, alpha = 0.01, delta = 0.01, burnin = 500, iter = 1000, thin = 100))
      tm <- topicmodels::posterior(mod)$topics
      full_tm <- matrix(NA, nrow = n_disease, ncol = k, dimnames = list(diseases, NULL))
      full_tm[rownames(tm), ] <- tm
      return(full_tm)
    }, error = function(e) return(NULL))
  }
  
  if (.Platform$OS.type == "unix") {
    res <- mclapply(1:n_iter, run_single, mc.cores = n_cores)
  } else {
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {library(topicmodels); library(tm)})
    clusterExport(cl, varlist=c("dtm_matrix","k","n_disease","diseases"), envir=environment())
    res <- parLapply(cl, 1:n_iter, run_single)
    stopCluster(cl)
  }
  return(res[!sapply(res, is.null)])
}

calculate_robust_similarity <- function(bootstrap_list, penalty_lambda = 1.0) {
  if(length(bootstrap_list) < 2) return(NULL)
  diseases <- rownames(bootstrap_list[[1]])
  n_boot <- length(bootstrap_list)
  sim_array <- array(NA, dim = c(length(diseases), length(diseases), n_boot))
  
  for (i in 1:n_boot) {
    tm <- bootstrap_list[[i]]; tm[is.na(tm)] <- 0
    sim <- cor(t(tm)); sim <- (sim + 1) / 2 
    sim_array[, , i] <- sim
  }
  
  mean_sim <- apply(sim_array, c(1, 2), mean, na.rm = TRUE)
  sd_sim <- apply(sim_array, c(1, 2), sd, na.rm = TRUE)
  robust_sim <- mean_sim - (penalty_lambda * sd_sim)
  robust_sim[robust_sim < 0] <- 0; diag(robust_sim) <- 1
  rownames(robust_sim) <- colnames(robust_sim) <- diseases
  
  disease_volatility <- rowMeans(sd_sim, na.rm = TRUE)
  conf <- 1 - (disease_volatility - min(disease_volatility))/(max(disease_volatility) - min(disease_volatility))
  names(conf) <- diseases
  return(list(robust_sim = robust_sim, confidence_score = conf))
}

############################################################################
# 5 & 6 & 7. Phase 1 Dispatch (Bootstrap OR Fast Mode)
############################################################################

if (is.null(opt$k_only)) {
  # =========================================================
  # 正常路线：执行 Bootstrap -> 选择最佳 K -> 画稳定性图
  # =========================================================
  results_by_k <- list() 
  stability_summary <- data.frame(k = integer(), mean_confidence = numeric(), threshold_used = numeric())
  
  cat("\n============================================\n")
  cat("Phase 1: Bootstrapping and Generating Stability Scores per K\n")
  cat("============================================\n")
  
  for (k in k_range) {
    cat(sprintf("\nProcessing K = %d ...\n", k))
    bs_results <- run_lda_bootstrap_parallel(input_mat, k = k, n_iter = n_bootstrap)
    
    if (length(bs_results) < 10) { cat("  [Warning] Too few successful runs. Skipping.\n"); next }
    res <- calculate_robust_similarity(bs_results, penalty_lambda = 1.0)
    
    if(!is.null(res)) {
      scores <- res$confidence_score
      avg_conf <- mean(scores, na.rm = TRUE)
      thresh_50 <- as.numeric(quantile(scores, 0.50, na.rm = TRUE))
      current_thresh <- min(0.25, thresh_50)
      
      stability_summary <- rbind(stability_summary, data.frame(k = as.numeric(k), mean_confidence = avg_conf, threshold_used = current_thresh))
      results_by_k[[as.character(k)]] <- res
      cat(sprintf("  -> K=%d | Mean Bootstrap Confidence: %.4f | Downstream Thresh: %.3f\n", k, avg_conf, current_thresh))
    }
  }
  saveRDS(list(summary = stability_summary, details = results_by_k), paste0(prefix, "_all_k_stability_results.rds"))
  
  if (nrow(stability_summary) == 0) stop("Stability summary is empty.")
  
  max_conf <- max(stability_summary$mean_confidence, na.rm = TRUE)
  tolerance <- 0.02 
  candidates <- stability_summary %>% filter(mean_confidence >= (max_conf - tolerance)) %>% arrange(k)
  
  best_k <- candidates$k[1]
  best_thresh <- candidates$threshold_used[1]
  cat(sprintf("\n>>> Best K Selected: %d <<<\n", best_k))
  
  k_key <- as.character(best_k)
  final_res <- results_by_k[[k_key]]
  if (is.null(names(final_res$confidence_score))) names(final_res$confidence_score) <- rownames(final_res$robust_sim)
  
  disease_ids <- rownames(final_res$robust_sim)
  meta_df <- data.frame(Disease = disease_ids, stringsAsFactors = FALSE)
  meta_df$clean_id <- gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$", "", meta_df$Disease)
  
  if (!is.null(master_map)) {
    master_sub <- master_map[!duplicated(master_map$clean_filename), ]
    meta_df <- merge(meta_df, master_sub, by.x = "clean_id", by.y = "clean_filename", all.x = TRUE)
    meta_df$Display_Name <- ifelse(!is.na(meta_df$trait) & meta_df$trait != "", meta_df$trait, meta_df$Disease)
    meta_df$label <- ifelse(!is.na(meta_df$label1) & meta_df$label1 != "", meta_df$label1, "Other")
    meta_df$label2 <- ifelse(!is.na(meta_df$label2) & meta_df$label2 != "", meta_df$label2, "Unknown")
  } else {
    meta_df$Display_Name <- meta_df$Disease; meta_df$label <- "Other"; meta_df$label2 <- "Unknown"
  }
  
  col_mean <- "#2E5C8A"; col_best <- "#E6550D"  
  p_combined <- ggplot(stability_summary, aes(x = k, y = mean_confidence)) +
    geom_line(color = col_mean, linewidth = 1.5) + geom_point(color = col_mean, size = 5) +
    geom_vline(xintercept = best_k, color = col_best, linetype = "dashed", linewidth = 1.5) +
    ggplot2::annotate("text", x = best_k + 0.15, y = min(stability_summary$mean_confidence, na.rm = TRUE), label = paste0("Best K = ", best_k), color = col_best, fontface = "bold", size = 6, hjust = 0, vjust = 0) +
    scale_x_continuous(breaks = k_range, labels = paste0(" ", k_range, " ")) + theme_classic(base_size = 18) +
    labs(title = "Optimal K Selection: Bootstrap Mean Confidence", x = "Number of Topics (K)", y = "Bootstrap Mean Confidence Score")
  ggsave(paste0(prefix, "_K_selection.png"), p_combined, width = 12, height = 7, dpi = 300)
  
  write.csv(stability_summary[stability_summary$k == best_k, ], paste0(prefix, "_best_k_selection_info.csv"), row.names = FALSE)
  
  scores <- final_res$confidence_score
  passed_ids <- if (opt$keep_all_traits) names(scores) else names(scores)[scores >= best_thresh]
  
  if(length(passed_ids) > 0) {
    stable_list_df <- meta_df[meta_df$Disease %in% passed_ids, ]
    stable_list_df$Confidence_Score <- scores[stable_list_df$Disease]
    cols_to_keep <- c("Disease", "Display_Name", "Confidence_Score", "label", "label2")
    stable_list_df <- stable_list_df[, cols_to_keep[cols_to_keep %in% colnames(stable_list_df)]]
    colnames(stable_list_df)[colnames(stable_list_df) == "label"] <- "Category" 
    colnames(stable_list_df)[colnames(stable_list_df) == "label2"] <- "Source"   
    write.csv(stable_list_df[order(-stable_list_df$Confidence_Score), ], paste0(prefix, "_stable_diseases_list.csv"), row.names = FALSE)
    
    plot_df <- stable_list_df
    plot_df$Display_Name_Unique <- make.unique(as.character(plot_df$Display_Name))
    plot_df <- plot_df[order(plot_df$Confidence_Score), ]
    plot_df$Display_Name_Unique <- factor(plot_df$Display_Name_Unique, levels = plot_df$Display_Name_Unique)
    
    custom_colors <- get_network_colors()
    unique_cats <- unique(plot_df$Category)
    cols <- setNames(rep(custom_colors, length.out=length(unique_cats)), unique_cats)
    shapes <- setNames(rep(c(16, 17, 15, 18, 3, 4, 8, 1, 2, 5), length.out=length(unique(plot_df$Source))), unique(plot_df$Source))
    plot_height <- max(6, nrow(plot_df) * 0.22)
    
    p_rank <- ggplot(plot_df, aes(x = Confidence_Score, y = Display_Name_Unique)) +
      geom_segment(aes(x = 0, xend = Confidence_Score, y = Display_Name_Unique, yend = Display_Name_Unique), color = "grey70", linewidth = 1) +
      geom_point(aes(color = Category, shape = Source), size = 4) +
      geom_vline(xintercept = best_thresh, color = "firebrick", linetype = "dashed", linewidth = 1) +
      ggplot2::annotate("text", x = best_thresh + 0.02, y = 1.5, label = paste0("Threshold (", round(best_thresh, 2), ")"), color = "firebrick", fontface = "bold", hjust = 0, size = 5) +
      scale_y_discrete(labels = setNames(as.character(plot_df$Display_Name), plot_df$Display_Name_Unique)) +
      scale_color_manual(values = cols) + scale_shape_manual(values = shapes) +
      scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
      theme_classic(base_size = 14) +
      theme(
        axis.text.y = element_text(face = "bold", size = 11, color = "black"), axis.text.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(), axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 15)),
        legend.position = "right", legend.title = element_text(face = "bold")
      ) +
      labs(title = paste0("Disease Stability (K = ", best_k, " )"), x = "Confidence Score")
    
    ggsave(paste0(prefix, "_ranking_best_k", best_k, ".jpg"), p_rank, width = 12, height = plot_height, dpi = 300, limitsize = FALSE)
  }
  
} else {
  # =========================================================
  # Fast Mode
  # =========================================================
  best_k <- opt$k_only
  opt$keep_all_traits <- TRUE # 隐式强制保留所有疾病
  best_thresh <- 0 
  
  cat("\n============================================\n")
  cat(sprintf("FAST MODE TRIGGERED: Skipping Bootstrapping. Forcing K = %d\n", best_k))
  cat("============================================\n")
  
  disease_ids <- rownames(input_mat)
  
  # 伪造全 1.0 的 Confidence Score 欺骗下游画图，防止报错
  dummy_conf <- rep(1.0, length(disease_ids))
  names(dummy_conf) <- disease_ids
  final_res <- list(confidence_score = dummy_conf)
  
  # 解析映射表以获取下游需要的 meta_df 格式
  meta_df <- data.frame(Disease = disease_ids, stringsAsFactors = FALSE)
  meta_df$clean_id <- gsub("^finngen_R12_|\\.list$|\\.tsv\\.bgz$", "", meta_df$Disease)
  
  if (!is.null(master_map)) {
    master_sub <- master_map[!duplicated(master_map$clean_filename), ]
    meta_df <- merge(meta_df, master_sub, by.x = "clean_id", by.y = "clean_filename", all.x = TRUE)
    meta_df$Display_Name <- ifelse(!is.na(meta_df$trait) & meta_df$trait != "", meta_df$trait, meta_df$Disease)
    meta_df$label <- ifelse(!is.na(meta_df$label1) & meta_df$label1 != "", meta_df$label1, "Other")
    meta_df$label2 <- ifelse(!is.na(meta_df$label2) & meta_df$label2 != "", meta_df$label2, "Unknown")
  } else {
    meta_df$Display_Name <- meta_df$Disease; meta_df$label <- "Other"; meta_df$label2 <- "Unknown"
  }
}

############################################################################
#        PART 2: Downstream Analysis Dispatch                              #
############################################################################

if (opt$mode == "gene") {
  if (opt$model %in% c("lda", "both")) run_downstream_lda(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, opt$seed, opt$keep_all_traits)
  if (opt$model %in% c("gaussian", "both")) run_downstream_gaussian(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, opt$seed, opt$keep_all_traits)
} else if (opt$mode == "ss") {
  if (opt$model %in% c("lda", "both")) run_downstream_ss_lda(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, opt$seed, opt$keep_all_traits)
  if (opt$model %in% c("gaussian", "both")) run_downstream_ss_gaussian(input_mat, final_res, meta_df, best_k, best_thresh, master_map, prefix, output_dir, opt$seed, opt$keep_all_traits)
}

cat("\nPipeline Execution Finished Successfully.\n")