# pathway_analysis.R
#' Pathway Enrichment Analysis for SNPic Models
#' 
#' This script performs comprehensive pathway enrichment analysis for both Normal and Gaussian SNPic models,
#' including GO Biological Process and KEGG pathway analysis with customizable gene numbers and visualization.
#' Updates: 
#' 1. Pagination max 8 topics/page (2 cols x 4 rows).
#' 2. Pathway names font size extremely enlarged for faceted plots.
#' 3. Facet strip headers show ONLY numbers (e.g., "1", "2") instead of "Topic 1".
#' 4. KEGG search drastically relaxed (Unadjusted P < 0.1, minGSSize=2, use_internal_data=FALSE)
#' 5. Faceted plots strict limit: Max 3 terms per topic, X-axis text hidden, Width expanded.
#' 6. Single plots: Max 5 terms by default, reduced to 3 if any text wraps.
#' 7. Added automatic CSV export for the top 3 pathways of each topic.

# Load required packages
load_pathway_packages <- function() {
  required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", 
                         "DOSE", "ggplot2", "dplyr", "tidyr", "openxlsx", "stringr")
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
    library(pkg, character.only = TRUE)
  }
  cat("All required packages loaded successfully.\n")
}

#' Extract top genes per topic from SNPic models
extract_top_genes_per_topic <- function(top_genes_obj, top_n = 50, model_type = "normal") {
  topic_genes <- list()
  
  if (model_type == "normal") {
    topics <- unique(top_genes_obj$Topic)
    # 按照数字顺序排序
    topics <- topics[order(as.numeric(topics))]
    
    for (topic in topics) {
      topic_data <- top_genes_obj %>% 
        filter(Topic == topic) %>%
        arrange(desc(Weight)) %>%
        head(top_n)
      
      topic_genes[[paste0("Topic", topic)]] <- topic_data$Comorbidity
    }
  } else if (model_type == "gaussian") {
    if (is.matrix(top_genes_obj)) {
      for (i in 1:ncol(top_genes_obj)) {
        topic_weights <- top_genes_obj[, i]
        top_indices <- order(topic_weights, decreasing = TRUE)[1:top_n]
        gene_names <- rownames(top_genes_obj)[top_indices]
        topic_genes[[paste0("Topic", i)]] <- gene_names
      }
    } else if (is.list(top_genes_obj)) {
      for (i in 1:length(top_genes_obj)) {
        topic_weights <- top_genes_obj[[i]]
        top_indices <- order(topic_weights, decreasing = TRUE)[1:top_n]
        gene_names <- names(topic_weights)[top_indices]
        topic_genes[[paste0("Topic", i)]] <- gene_names
      }
    }
  }
  
  cat("Extracted top", top_n, "genes for", length(topic_genes), "topics from", model_type, "model\n")
  return(topic_genes)
}

#' Convert gene symbols to Entrez IDs
symbol_to_entrez <- function(gene_symbols) {
  clean_genes <- gsub("\\..*$", "", gene_symbols)
  clean_genes <- gsub("_", "-", clean_genes)
  
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = clean_genes,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  return(entrez_ids)
}

#' Perform enrichment analysis for a gene list
perform_enrichment_analysis <- function(gene_list, analysis_name, p_cutoff = 0.05, q_cutoff = 0.2, min_genes = 5) {
  cat("Enrichment analysis for", analysis_name, "with", length(gene_list), "genes\n")
  
  entrez_ids <- symbol_to_entrez(gene_list)
  entrez_ids <- na.omit(entrez_ids)
  
  if (length(entrez_ids) < min_genes) {
    cat("Warning: Only", length(entrez_ids), "valid ENTREZ IDs found - skipping\n")
    return(NULL)
  }
  
  results <- list()
  
  # 极度放宽 KEGG 搜索能力 (使用原始 P 值 < 0.1)
  tryCatch({
    kegg_result <- enrichKEGG(
      gene = entrez_ids,
      organism = 'hsa',
      pAdjustMethod = "none",  
      pvalueCutoff = 0.1,      
      qvalueCutoff = 1.0,      
      minGSSize = 2,           
      use_internal_data = FALSE 
    )
    results$kegg <- kegg_result
    cat("  - KEGG:", ifelse(!is.null(kegg_result), nrow(kegg_result), 0), "pathways found (Unadjusted P < 0.1)\n")
  }, error = function(e) {
    cat("  - KEGG analysis failed:", e$message, "\n")
    results$kegg <- NULL
  })
  
  # GO Biological Process enrichment (保持原有标准)
  tryCatch({
    go_bp_result <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = p_cutoff,
      pAdjustMethod = "BH",
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
    results$go_bp <- go_bp_result
    cat("  - GO BP:", ifelse(!is.null(go_bp_result), nrow(go_bp_result), 0), "terms found\n")
  }, error = function(e) {
    cat("  - GO BP analysis failed:", e$message, "\n")
    results$go_bp <- NULL
  })
  
  return(results)
}

#' Generate enrichment plots for a topic (Single topic plot)
generate_enrichment_plots <- function(enrichment_results, model_name, topic_name, prefix_pathway) {
  plot_dir <- file.path(prefix_pathway, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  text_width_single <- 40 
  
  # KEGG plots
  if (!is.null(enrichment_results$kegg) && nrow(enrichment_results$kegg) > 0) {
    tryCatch({
      temp_df <- head(as.data.frame(enrichment_results$kegg), 5)
      temp_df$Description <- str_wrap(temp_df$Description, width = text_width_single)
      
      show_cat <- if (any(str_detect(temp_df$Description, "\n"))) 3 else 5
      
      enrichment_results$kegg@result$Description <- str_wrap(enrichment_results$kegg@result$Description, width = text_width_single)
      
      p1 <- dotplot(enrichment_results$kegg, showCategory = show_cat, color = "pvalue") +
        ggtitle(paste("KEGG Enrichment -", model_name, topic_name)) +
        theme_minimal(base_size = 20) +
        theme(
          axis.text.y = element_text(size = 24, face = "bold", color = "black", lineheight = 0.8),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 22, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20))
        ) 
      
      ggsave(file.path(plot_dir, paste0(model_name, "_", topic_name, "_kegg_dotplot.png")), 
             p1, width = 14, height = 10, dpi = 300, bg = "white")
    }, error = function(e) NULL)
  }
  
  # GO BP plots
  if (!is.null(enrichment_results$go_bp) && nrow(enrichment_results$go_bp) > 0) {
    tryCatch({
      temp_df <- head(as.data.frame(enrichment_results$go_bp), 5)
      temp_df$Description <- str_wrap(temp_df$Description, width = text_width_single)
      
      show_cat <- if (any(str_detect(temp_df$Description, "\n"))) 3 else 5
      
      enrichment_results$go_bp@result$Description <- str_wrap(enrichment_results$go_bp@result$Description, width = text_width_single)
      
      p2 <- dotplot(enrichment_results$go_bp, showCategory = show_cat) +
        ggtitle(paste("GO Biological Process -", model_name, topic_name)) +
        theme_minimal(base_size = 20) +
        theme(
          axis.text.y = element_text(size = 24, face = "bold", color = "black", lineheight = 0.8),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.title = element_text(size = 22, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 20))
        )
        
      ggsave(file.path(plot_dir, paste0(model_name, "_", topic_name, "_go_bp_dotplot.png")), 
             p2, width = 14, height = 10, dpi = 300, bg = "white")
    }, error = function(e) NULL)
  }
}

#' Analyze all topics for a model
analyze_all_topics <- function(topic_genes_list, model_name, prefix_pathway, top_n = 50) {
  all_results <- list()
  
  topic_names <- names(topic_genes_list)
  topic_nums <- as.numeric(gsub("Topic", "", topic_names))
  sorted_indices <- order(topic_nums)
  sorted_topic_names <- topic_names[sorted_indices]
  
  for (topic_name in sorted_topic_names) {
    cat("Analyzing", topic_name, "from", model_name, "model\n")
    
    genes <- topic_genes_list[[topic_name]]
    if (length(genes) > top_n) {
      genes <- genes[1:top_n]  
    }
    
    enrichment_results <- perform_enrichment_analysis(genes, paste0(model_name, "_", topic_name))
    
    if (!is.null(enrichment_results)) {
      all_results[[topic_name]] <- enrichment_results
      generate_enrichment_plots(enrichment_results, model_name, topic_name, prefix_pathway)
    } else {
      all_results[[topic_name]] <- list(kegg = NULL, go_bp = NULL)
    }
  }
  
  return(all_results)
}

#' 【新增功能】：将每个 Topic 最显著的 Pathways 导出成表格
#' @param enrichment_results Enrichment results list
#' @param model_name Name of the model
#' @param prefix_pathway Output directory path
#' @param top_n Number of top terms to extract per topic
export_top_pathways_table <- function(enrichment_results, model_name, prefix_pathway, top_n = 3) {
  all_data <- data.frame()
  
  topic_names <- names(enrichment_results)
  topic_nums <- as.numeric(gsub("Topic", "", topic_names))
  sorted_topics <- topic_names[order(topic_nums)]
  
  for (topic in sorted_topics) {
    # 提取 GO BP 前三名
    if (!is.null(enrichment_results[[topic]]$go_bp) && nrow(as.data.frame(enrichment_results[[topic]]$go_bp)) > 0) {
      go_df <- as.data.frame(enrichment_results[[topic]]$go_bp)
      go_df <- go_df %>%
        arrange(p.adjust) %>%
        head(top_n) %>%
        mutate(Topic = topic, Model = model_name, Analysis = "GO_BP") %>%
        dplyr::select(Topic, Model, Analysis, ID, Description, pvalue, p.adjust, Count)
      all_data <- rbind(all_data, go_df)
    } else {
      # 若无显著结果则填入占位符
      all_data <- rbind(all_data, data.frame(Topic = topic, Model = model_name, Analysis = "GO_BP", 
                                             ID = NA, Description = "No significant enrichment", 
                                             pvalue = NA, p.adjust = NA, Count = 0))
    }
    
    # 提取 KEGG 前三名
    if (!is.null(enrichment_results[[topic]]$kegg) && nrow(as.data.frame(enrichment_results[[topic]]$kegg)) > 0) {
      kegg_df <- as.data.frame(enrichment_results[[topic]]$kegg)
      kegg_df <- kegg_df %>%
        arrange(pvalue) %>% # KEGG 使用未校正的 pvalue 排序
        head(top_n) %>%
        mutate(Topic = topic, Model = model_name, Analysis = "KEGG") %>%
        dplyr::select(Topic, Model, Analysis, ID, Description, pvalue, p.adjust, Count)
      all_data <- rbind(all_data, kegg_df)
    } else {
      all_data <- rbind(all_data, data.frame(Topic = topic, Model = model_name, Analysis = "KEGG", 
                                             ID = NA, Description = "No significant enrichment", 
                                             pvalue = NA, p.adjust = NA, Count = 0))
    }
  }
  
  if (nrow(all_data) > 0) {
    out_file <- file.path(prefix_pathway, paste0(tolower(model_name), "_top", top_n, "_pathways_summary.csv"))
    write.csv(all_data, out_file, row.names = FALSE)
    cat(sprintf("  - Successfully exported top %d pathways summary table to: %s\n", top_n, out_file))
  }
}


#' Create combined GO BP and KEGG faceted dotplots with Pagination (2 cols x 4 rows)
#' 
#' @param enrichment_results Enrichment results
#' @param model_name Name of the model
#' @param prefix_pathway Output directory path
#' @param top_terms Number of top terms to show per topic (default: 3)
create_combined_faceted_plots <- function(enrichment_results, model_name, prefix_pathway, top_terms = 3) {
  
  # Function to extract and combine enrichment data
  extract_enrichment_data <- function(enrichment_results, analysis_type) {
    all_data <- data.frame()
    
    topic_names <- names(enrichment_results)
    topic_nums <- as.numeric(gsub("Topic", "", topic_names))
    sorted_indices <- order(topic_nums)
    sorted_topic_names <- topic_names[sorted_indices]
    
    sample_structure <- NULL
    for (topic_name in sorted_topic_names) {
      if (!is.null(enrichment_results[[topic_name]][[analysis_type]]) && 
          nrow(enrichment_results[[topic_name]][[analysis_type]]) > 0) {
        sample_structure <- as.data.frame(enrichment_results[[topic_name]][[analysis_type]])
        break
      }
    }
    
    for (topic_name in sorted_topic_names) {
      if (!is.null(enrichment_results[[topic_name]][[analysis_type]]) && 
          nrow(enrichment_results[[topic_name]][[analysis_type]]) > 0) {
        
        enrich_df <- as.data.frame(enrichment_results[[topic_name]][[analysis_type]])
        enrich_df$Topic <- topic_name
        enrich_df$Model <- model_name
        enrich_df$Analysis <- analysis_type
        
        # 初始提取，强制只取 top_terms (现在被设定为 3)
        if (analysis_type == "kegg") {
           top_data <- enrich_df %>%
             arrange(pvalue) %>%
             head(top_terms)
        } else {
           top_data <- enrich_df %>%
             arrange(p.adjust) %>%
             head(top_terms)
        }
        
        all_data <- rbind(all_data, top_data)
      } else {
        # 处理空结果
        if (!is.null(sample_structure)) {
          empty_row <- sample_structure[1, , drop = FALSE]
          empty_row[1, ] <- NA
          empty_row$Topic <- topic_name
          empty_row$Model <- model_name
          empty_row$Analysis <- analysis_type
          empty_row$Description <- "No significant enrichment"
          empty_row$Count <- 0
          empty_row$pvalue <- 1.0
          empty_row$p.adjust <- 1.0
          all_data <- rbind(all_data, empty_row)
        } else {
          empty_row <- data.frame(
            ID = NA, Description = "No significant enrichment",
            GeneRatio = NA, BgRatio = NA, pvalue = 1.0, p.adjust = 1.0, 
            qvalue = NA, geneID = NA, Count = 0,
            Topic = topic_name, Model = model_name, Analysis = analysis_type,
            stringsAsFactors = FALSE
          )
          all_data <- rbind(all_data, empty_row)
        }
      }
    }
    return(all_data)
  }
  
  combined_plot_dir <- file.path(prefix_pathway, "combined_plots")
  dir.create(combined_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Function to create faceted dotplot
  create_faceted_dotplot <- function(plot_data, analysis_name, page_suffix = "") {
    if (nrow(plot_data) == 0 || all(plot_data$Description == "No significant enrichment")) {
      cat("No significant", analysis_name, "enrichment available for", model_name, page_suffix, "\n")
      return(NULL)
    }
    
    valid_topics <- plot_data %>%
      group_by(Topic) %>%
      filter(!all(Description == "No significant enrichment")) %>%
      pull(Topic) %>%
      unique()
    
    num_valid_topics <- length(valid_topics)
    if (num_valid_topics == 0) return(NULL)
    
    current_ncol <- 2  
    if (num_valid_topics < 2) current_ncol <- 1
    
    text_width <- 60 
    dot_range <- c(4, 10)
    
    plot_data_processed <- plot_data %>%
      filter(Topic %in% valid_topics) %>%
      mutate(
        Topic_Num = as.numeric(gsub("Topic", "", Topic)),
        Topic_Label = as.character(Topic_Num),
        DescriptionRaw = ifelse(is.na(Description), "No significant enrichment", Description),
        DescriptionWrapped = str_wrap(DescriptionRaw, width = text_width),
        HasWrap = str_detect(DescriptionWrapped, "\n"),
        Topic = factor(Topic_Label, levels = unique(as.character(sort(unique(Topic_Num))))),
        
        # 计算对数P值用于X轴
        neg_log_pvalue = if(analysis_name == "kegg") {
                            ifelse(is.na(pvalue) | pvalue == 1.0, 0, -log10(pvalue))
                         } else {
                            ifelse(is.na(p.adjust) | p.adjust == 1.0, 0, -log10(p.adjust))
                         },
        Count = ifelse(is.na(Count), 0, Count)
      ) %>%
      group_by(Topic) %>%
      arrange(desc(neg_log_pvalue)) %>%
      mutate(
        CurrentSliceLimit = top_terms, # 【修改点】：动态跟随顶层参数，不再硬编码为 2
        GroupRank = row_number()
      ) %>%
      filter(GroupRank <= CurrentSliceLimit) %>%
      mutate(Description = DescriptionWrapped) %>%
      dplyr::select(-DescriptionRaw, -DescriptionWrapped, -HasWrap, -CurrentSliceLimit, -GroupRank) %>%
      ungroup()
    
    base_font_size <- 24 
    
    x_label <- ifelse(analysis_name == "kegg", "-log10(pvalue)", "-log10(p.adjust)")
    
    p <- ggplot(plot_data_processed, aes(x = neg_log_pvalue, y = reorder(Description, neg_log_pvalue))) +
      geom_point(aes(size = Count, color = neg_log_pvalue)) +
      
      facet_wrap(~ Topic, scales = "free", ncol = current_ncol) + 
      
      scale_color_gradient(low = "blue", high = "red", name = x_label) +
      scale_size_continuous(name = "Count", range = dot_range) + 
      
      theme_minimal(base_size = base_font_size) + 
      theme(
        axis.text.y = element_text(size = 36, face = "bold", lineheight = 0.8, color = "black"), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(size = 28, face = "bold", margin = margin(t = 10)),   
        axis.title.y = element_blank(),
        
        strip.text = element_text(face = "bold", size = 32, margin = margin(12, 0, 12, 0)), 
        strip.background = element_rect(fill = "grey92", color = NA),
        
        legend.title = element_text(size = 26, face = "bold"), 
        legend.text = element_text(size = 24),                 
        legend.position = "right",
        
        plot.title = element_text(hjust = 0.5, size = 40, face = "bold", margin = margin(15, 0, 30, 0)), 
        
        panel.spacing.x = unit(2.5, "lines"), 
        panel.spacing.y = unit(2.0, "lines"), 
        plot.margin = margin(40, 40, 40, 40),
        panel.grid.major.y = element_line(size = 0.3, color = "grey85", linetype = "dashed"), 
        panel.grid.major.x = element_blank() 
      ) +
      labs(
        title = paste(toupper(analysis_name), "Enrichment -", model_name, page_suffix),
        x = x_label,
        y = NULL
      )
    
    return(p)
  }
  
  save_plot <- function(p, filename, n_topics) {
    w <- 45
    
    # 【修改点】：高度统一增加 2
    if(n_topics > 6) { h <- 26 } 
    else if (n_topics > 4) { h <- 20 } 
    else if (n_topics > 2) { h <- 15 } 
    else { h <- 10 }
    
    ggsave(file.path(combined_plot_dir, filename),
           p, width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  }
  
  process_paginated_plots <- function(full_data, analysis_type) {
    if (nrow(full_data) == 0) return()
    
    all_topics <- unique(full_data$Topic)
    topic_nums <- as.numeric(gsub("Topic", "", all_topics))
    sorted_topics <- all_topics[order(topic_nums)]
    
    chunk_size <- 8
    total_topics <- length(sorted_topics)
    total_pages <- ceiling(total_topics / chunk_size)
    
    cat(sprintf("Processing %s: %d topics total, split into %d plots (max 8/plot)\n", 
                analysis_type, total_topics, total_pages))
    
    for (i in 1:total_pages) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, total_topics)
      
      current_topics <- sorted_topics[start_idx:end_idx]
      current_data <- full_data %>% filter(Topic %in% current_topics)
      
      if (total_pages > 1) {
        page_suffix <- paste0(" (Part ", i, "/", total_pages, ")")
        file_suffix <- paste0("_part", i)
      } else {
        page_suffix <- ""
        file_suffix <- ""
      }
      
      p <- create_faceted_dotplot(current_data, analysis_type, page_suffix)
      
      if (!is.null(p)) {
        n_current_topics <- length(unique(current_data$Topic))
        filename <- paste0(tolower(model_name), "_faceted_", analysis_type, "_dotplot", file_suffix, ".png")
        save_plot(p, filename, n_current_topics)
        cat(sprintf("  - Saved %s\n", filename))
      }
    }
  }
  
  # Process GO BP with pagination
  go_bp_data <- extract_enrichment_data(enrichment_results, "go_bp")
  process_paginated_plots(go_bp_data, "go_bp")
  
  # Process KEGG with pagination
  kegg_data <- extract_enrichment_data(enrichment_results, "kegg")
  process_paginated_plots(kegg_data, "kegg")
}

#' Main function to run pathway enrichment analysis
run_pathway_analysis <- function(normal_top_genes, gaussian_top_genes,
                                 normal_topic_distribution, gaussian_topic_distribution,
                                 prefix_pathway, top_genes = 50, top_terms = 3) { # 【修改点】：默认提取 3 个 terms
  
  load_pathway_packages()
  dir.create(prefix_pathway, recursive = TRUE, showWarnings = FALSE)
  
  cat("Starting pathway enrichment analysis...\n")
  cat("Output directory:", prefix_pathway, "\n")
  
  # Extract top genes
  normal_topic_genes <- extract_top_genes_per_topic(normal_top_genes, top_n = top_genes, model_type = "normal")
  gaussian_topic_genes <- extract_top_genes_per_topic(gaussian_top_genes, top_n = top_genes, model_type = "gaussian")
  
  # Save gene lists
  for (topic_name in names(normal_topic_genes)) {
    genes <- as.character(normal_topic_genes[[topic_name]])
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) > 0) {
      write.table(genes, file = paste0(prefix_pathway, "/normal_", topic_name, "_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  for (topic_name in names(gaussian_topic_genes)) {
    genes <- as.character(gaussian_topic_genes[[topic_name]])
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) > 0) {
      write.table(genes, file = paste0(prefix_pathway, "/gaussian_", topic_name, "_genes.txt"),
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  cat("Gene lists saved to:", prefix_pathway, "\n")
  
  # Perform enrichment analysis
  cat("Starting enrichment analysis for Normal SNPic model...\n")
  normal_enrichment_results <- analyze_all_topics(normal_topic_genes, "normal", prefix_pathway, top_n = top_genes)
  
  cat("Starting enrichment analysis for Gaussian SNPic model...\n")
  gaussian_enrichment_results <- analyze_all_topics(gaussian_topic_genes, "gaussian", prefix_pathway, top_n = top_genes)
  
  cat("Exporting top 3 pathways summary tables...\n")
  export_top_pathways_table(normal_enrichment_results, "Normal_SNPic", prefix_pathway, top_n = 3)
  export_top_pathways_table(gaussian_enrichment_results, "Gaussian_SNPic", prefix_pathway, top_n = 3)
  
  # Create combined faceted plots
  cat("Creating combined faceted plots (paginated max 8 topics/plot)...\n")
  create_combined_faceted_plots(normal_enrichment_results, "Normal_SNPic", prefix_pathway, top_terms = 3) # 【修改点】：调用时显式传入 3
  create_combined_faceted_plots(gaussian_enrichment_results, "Gaussian_SNPic", prefix_pathway, top_terms = 3) # 【修改点】：调用时显式传入 3
  
  # Save all results
  save(normal_topic_genes, gaussian_topic_genes,
       normal_enrichment_results, gaussian_enrichment_results,
       file = file.path(prefix_pathway, "all_enrichment_results.RData"))
  
  results <- list(
    normal_topic_genes = normal_topic_genes,
    gaussian_topic_genes = gaussian_topic_genes,
    normal_enrichment_results = normal_enrichment_results,
    gaussian_enrichment_results = gaussian_enrichment_results
  )
  
  cat("Pathway enrichment analysis completed successfully!\n")
  return(results)
}