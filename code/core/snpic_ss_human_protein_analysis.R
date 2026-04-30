#' Human Gene Tissue Expression Analysis (GTEx v8 Full 54 Tissues - REAL DATA ONLY)
#' 
#' Connect to GTEx database, map Gene Symbols to ENSEMBL IDs, fetch real tissue 
#' expression data across 54 tissues, and generate comprehensive visualizations.
#' 

library(dplyr)
library(ggplot2)
library(plotly)
library(httr)
library(jsonlite)
library(networkD3)
library(tidyr)
library(stringr)
library(ggalluvial)
library(RColorBrewer)
library(viridis)
library(scales)
library(tidytext) 

# GTEx API configuration
GTEX_BASE_URL <- "https://gtexportal.org/api/v2/"
GTEX_TIMEOUT <- 30
GTEX_RETRY_ATTEMPTS <- 3

# GTEx API v2 to standard display name mapping
gtex_tissue_mapping <- c(
  "Brain_Amygdala" = "Brain - Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24" = "Brain - Anterior cingulate cortex (BA24)",
  "Brain_Caudate_basal_ganglia" = "Brain - Caudate (basal ganglia)",
  "Brain_Cerebellar_Hemisphere" = "Brain - Cerebellar Hemisphere",
  "Brain_Cerebellum" = "Brain - Cerebellum",
  "Brain_Cortex" = "Brain - Cortex",
  "Brain_Frontal_Cortex_BA9" = "Brain - Frontal Cortex (BA9)",
  "Brain_Hippocampus" = "Brain - Hippocampus",
  "Brain_Hypothalamus" = "Brain - Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia" = "Brain - Nucleus accumbens (basal ganglia)",
  "Brain_Putamen_basal_ganglia" = "Brain - Putamen (basal ganglia)",
  "Brain_Spinal_cord_cervical_c-1" = "Brain - Spinal cord (cervical c-1)",
  "Brain_Substantia_nigra" = "Brain - Substantia nigra",
  "Nerve_Tibial" = "Nerve - Tibial",
  "Heart_Atrial_Appendage" = "Heart - Atrial Appendage",
  "Heart_Left_Ventricle" = "Heart - Left Ventricle",
  "Artery_Aorta" = "Artery - Aorta",
  "Artery_Coronary" = "Artery - Coronary",
  "Artery_Tibial" = "Artery - Tibial",
  "Esophagus_Gastroesophageal_Junction" = "Esophagus - Gastroesophageal Junction",
  "Esophagus_Mucosa" = "Esophagus - Mucosa",
  "Esophagus_Muscularis" = "Esophagus - Muscularis",
  "Stomach" = "Stomach",
  "Small_Intestine_Terminal_Ileum" = "Small Intestine - Terminal Ileum",
  "Colon_Sigmoid" = "Colon - Sigmoid",
  "Colon_Transverse" = "Colon - Transverse",
  "Liver" = "Liver",
  "Pancreas" = "Pancreas",
  "Thyroid" = "Thyroid",
  "Adrenal_Gland" = "Adrenal Gland",
  "Pituitary" = "Pituitary",
  "Kidney_Cortex" = "Kidney - Cortex",
  "Kidney_Medulla" = "Kidney - Medulla",
  "Whole_Blood" = "Whole Blood",
  "Spleen" = "Spleen",
  "Cells_EBV-transformed_lymphocytes" = "Cells - EBV-transformed lymphocytes",
  "Skin_Not_Sun_Exposed_Suprapubic" = "Skin - Not Sun Exposed (Suprapubic)",
  "Skin_Sun_Exposed_Lower_leg" = "Skin - Sun Exposed (Lower leg)",
  "Adipose_Subcutaneous" = "Adipose - Subcutaneous",
  "Adipose_Visceral_Omentum" = "Adipose - Visceral (Omentum)",
  "Cells_Cultured_fibroblasts" = "Cells - Cultured fibroblasts",
  "Breast_Mammary_Tissue" = "Breast - Mammary Tissue",
  "Ovary" = "Ovary",
  "Uterus" = "Uterus",
  "Vagina" = "Vagina",
  "Fallopian_Tube" = "Fallopian Tube",
  "Cervix_Ectocervix" = "Cervix - Ectocervix",
  "Cervix_Endocervix" = "Cervix - Endocervix",
  "Prostate" = "Prostate",
  "Testis" = "Testis",
  "Muscle_Skeletal" = "Muscle - Skeletal",
  "Lung" = "Lung",
  "Bladder" = "Bladder",
  "Minor_Salivary_Gland" = "Minor Salivary Gland"
)

#' Clean a single gene name
clean_single_gene <- function(gene) {
  if (is.na(gene) || is.null(gene)) return(NA)
  gene <- str_replace(gene, "\\(.*\\)", "")
  if (str_detect(gene, ",")) gene <- str_split(gene, ",")[[1]][1]
  return(str_trim(gene))
}

clean_gene_names <- function(gene_names) sapply(gene_names, clean_single_gene, USE.NAMES = FALSE)

extract_unique_genes <- function(gene_names) {
  cleaned_genes <- clean_gene_names(gene_names)
  unique(cleaned_genes[!is.na(cleaned_genes)])
}

#' Validate gene symbols and fetch their ENSEMBL (gencodeId) mappings from GTEx
validate_and_map_genes <- function(gene_symbols) {
  cat("Validating gene symbols and fetching ENSEMBL IDs from GTEx...\n")
  endpoint <- "reference/gene"
  url <- paste0(GTEX_BASE_URL, endpoint)
  
  valid_mapping <- list()
  
  for (gene in gene_symbols) {
    query_params <- list(format = "json", geneId = gene)
    
    tryCatch({
      response <- httr::GET(url, query = query_params, httr::timeout(10))
      if (httr::status_code(response) == 200) {
        content_text <- httr::content(response, as = "text", encoding = "UTF-8")
        data <- jsonlite::fromJSON(content_text)
        
        if (!is.null(data$data) && nrow(data$data) > 0) {
          match_idx <- which(toupper(data$data$geneSymbol) == toupper(gene))
          if (length(match_idx) > 0) {
            valid_mapping[[gene]] <- data$data$gencodeId[match_idx[1]]
          }
        }
      }
    }, error = function(e) {})
    Sys.sleep(0.1) 
  }
  
  cat(sprintf("Successfully mapped %d out of %d genes to ENSEMBL IDs.\n", length(valid_mapping), length(gene_symbols)))
  return(valid_mapping)
}

#' Get gene expression data from GTEx portal using ENSEMBL IDs
get_gtex_expression <- function(gene_mapping) {
  cat("Fetching median expression data from GTEx database (This might take 1-2 minutes)...\n")
  
  all_expression_data <- list()
  genes <- names(gene_mapping)
  
  for (i in seq_along(genes)) {
    gene_symbol <- genes[i]
    gencode_id <- gene_mapping[[gene_symbol]]
    
    cat(sprintf("  Fetching %d/%d: %s (%s)...\n", i, length(genes), gene_symbol, gencode_id))
    
    endpoint <- "expression/medianGeneExpression"
    url <- paste0(GTEX_BASE_URL, endpoint)
    
    query_params <- list(
      format = "json",
      datasetId = "gtex_v8",
      gencodeId = gencode_id
    )
    
    for (attempt in 1:GTEX_RETRY_ATTEMPTS) {
      tryCatch({
        response <- httr::GET(url, query = query_params, httr::timeout(GTEX_TIMEOUT))
        
        if (httr::status_code(response) == 200) {
          content_text <- httr::content(response, as = "text", encoding = "UTF-8")
          data <- jsonlite::fromJSON(content_text)
          
          if (!is.null(data$data) && nrow(data$data) > 0) {
            expression_df <- data$data
            expression_df$geneSymbol <- gene_symbol 
            
            # Map tissueSiteDetailId to the expected tissueSiteDetail format
            if ("tissueSiteDetailId" %in% colnames(expression_df)) {
              expression_df$tissueSiteDetail <- gtex_tissue_mapping[expression_df$tissueSiteDetailId]
              expression_df$tissueSiteDetail[is.na(expression_df$tissueSiteDetail)] <- expression_df$tissueSiteDetailId[is.na(expression_df$tissueSiteDetail)]
            }
            
            all_expression_data[[gene_symbol]] <- expression_df
            break  
          }
        }
      }, error = function(e) {})
      if (attempt < GTEX_RETRY_ATTEMPTS) Sys.sleep(1)
    }
    Sys.sleep(0.1)
  }
  
  if (length(all_expression_data) > 0) {
    return(bind_rows(all_expression_data))
  } else {
    return(NULL)
  }
}

#' Prepare data for visualization with custom full 54 tissue ordering
prepare_visualization_data <- function(topic_tissue_data) {
  tissue_categories <- list(
    "Brain & Nerve" = c("Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)", "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra", "Nerve - Tibial"),
    "Cardiovascular" = c("Heart - Atrial Appendage", "Heart - Left Ventricle", "Artery - Aorta", "Artery - Coronary", "Artery - Tibial"),
    "Digestive" = c("Esophagus - Gastroesophageal Junction", "Esophagus - Mucosa", "Esophagus - Muscularis", "Stomach", "Small Intestine - Terminal Ileum", "Colon - Sigmoid", "Colon - Transverse"),
    "Endocrine & Metabolic" = c("Liver", "Pancreas", "Thyroid", "Adrenal Gland", "Pituitary", "Kidney - Cortex", "Kidney - Medulla"),
    "Immune & Blood" = c("Whole Blood", "Spleen", "Cells - EBV-transformed lymphocytes"),
    "Skin & Adipose" = c("Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)", "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Cells - Cultured fibroblasts"),
    "Reproductive" = c("Breast - Mammary Tissue", "Ovary", "Uterus", "Vagina", "Fallopian Tube", "Cervix - Ectocervix", "Cervix - Endocervix", "Prostate", "Testis"),
    "Musculoskeletal" = c("Muscle - Skeletal"),
    "Respiratory" = c("Lung"),
    "Other" = c("Bladder", "Minor Salivary Gland")
  )
  
  tissue_order <- unlist(tissue_categories)
  
  processed_data <- topic_tissue_data$full_data %>%
    mutate(
      Topic = paste("Topic", as.character(Topic)),
      tissueSiteDetail = factor(tissueSiteDetail, levels = tissue_order),
      tissueCategory = case_when(
        tissueSiteDetail %in% tissue_categories$`Brain & Nerve` ~ "Brain & Nerve",
        tissueSiteDetail %in% tissue_categories$Cardiovascular ~ "Cardiovascular",
        tissueSiteDetail %in% tissue_categories$Digestive ~ "Digestive",
        tissueSiteDetail %in% tissue_categories$`Endocrine & Metabolic` ~ "Endocrine & Metabolic",
        tissueSiteDetail %in% tissue_categories$`Immune & Blood` ~ "Immune & Blood",
        tissueSiteDetail %in% tissue_categories$`Skin & Adipose` ~ "Skin & Adipose",
        tissueSiteDetail %in% tissue_categories$Reproductive ~ "Reproductive",
        tissueSiteDetail %in% tissue_categories$Musculoskeletal ~ "Musculoskeletal",
        tissueSiteDetail %in% tissue_categories$Respiratory ~ "Respiratory",
        TRUE ~ "Other"
      )
    ) %>%
    arrange(tissueCategory, tissueSiteDetail)
  
  return(processed_data)
}

#' [核心修改区]：Analyze topic-tissue associations using Tissue Specificity Normalization
analyze_topic_tissue_association <- function(topic_gene_df, expression_data, top_tissues_per_topic = 5) {
  
  topic_gene_clean <- topic_gene_df %>% mutate(geneSymbol = clean_gene_names(Comorbidity))
  
  # 1. 基因层面的归一化：计算基因表达在各个组织的分配比例 (Gene-wise Specificity)
  expression_normalized <- expression_data %>%
    group_by(geneSymbol) %>%
    mutate(
      total_expr = sum(median, na.rm = TRUE),
      specificity_score = ifelse(total_expr > 0, median / total_expr, 0)
    ) %>%
    ungroup()
    
  # 2. [新增] 组织层面的基线背景计算 (Tissue-wise Background)
  # 计算在所有被查询的基因中，每个器官的平均背景表达比例（用于抓出睾丸这类广泛表达的刺客）
  tissue_background <- expression_normalized %>%
    group_by(tissueSiteDetail) %>%
    summarise(
      tissue_bg_score = mean(specificity_score, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # 3. 汇总至 Topic 层级，并进行背景相对倍数折算 (Fold Enrichment)
  topic_tissue_expression <- topic_gene_clean %>%
    left_join(expression_normalized, by = "geneSymbol", relationship = "many-to-many") %>%
    group_by(Topic, tissueSiteDetail) %>%
    summarise(
      mean_expression = mean(specificity_score, na.rm = TRUE), 
      raw_tpm_mean = mean(median, na.rm = TRUE), 
      n_genes = n_distinct(geneSymbol),
      .groups = 'drop'
    ) %>%
    filter(!is.na(mean_expression)) %>%
    # 联结组织背景分数
    left_join(tissue_background, by = "tissueSiteDetail") %>%
    mutate(
      # [核心算法] 相对富集分 = 该Topic平均分 / 组织的背景均分 (Fold Change)
      normalized_score = ifelse(tissue_bg_score > 0, mean_expression / tissue_bg_score, 0)
    ) %>%
    arrange(Topic, desc(normalized_score))
  
  top_tissues <- topic_tissue_expression %>%
    group_by(Topic) %>%
    # 改为依据归一化后的分数进行选取和排行
    slice_max(order_by = normalized_score, n = top_tissues_per_topic, with_ties = FALSE) %>%
    ungroup()
  
  return(list(full_data = topic_tissue_expression, top_tissues = top_tissues))
}

prepare_sankey_data <- function(topic_tissue_data, min_expression = 0) {
  topics_ordered <- sort(unique(as.numeric(as.character(topic_tissue_data$top_tissues$Topic))))
  tissues <- unique(topic_tissue_data$top_tissues$tissueSiteDetail)
  
  topic_nodes <- paste("Topic", topics_ordered)
  nodes <- data.frame(
    name = c(topic_nodes, tissues),
    type = c(rep("topic", length(topic_nodes)), rep("tissue", length(tissues)))
  )
  
  links <- topic_tissue_data$top_tissues %>%
    filter(normalized_score >= min_expression) %>%
    mutate(
      source = match(as.numeric(as.character(Topic)), topics_ordered) - 1,
      target = match(tissueSiteDetail, tissues) + length(topic_nodes) - 1,
      value = normalized_score # 更新为使用新的富集度
    ) %>%
    dplyr::select(source, target, value)
  
  return(list(links = links, nodes = nodes))
}

create_sankey_plot <- function(sankey_data, width = 1200, height = 800) {
  sankeyNetwork(
    Links = sankey_data$links, Nodes = sankey_data$nodes,
    Source = "source", Target = "target", Value = "value", NodeID = "name",
    units = "Enrichment", fontSize = 16, nodeWidth = 30, nodePadding = 15,
    width = width, height = height, NodeGroup = "type",
    colourScale = JS('d3.scaleOrdinal().domain(["topic", "tissue"]).range(["#1f77b4", "#ff7f0e"])')
  )
}

create_static_sankey_plot <- function(topic_tissue_data, color_palette = "Set3", alpha = 0.7) {
  sankey_data <- topic_tissue_data$top_tissues %>%
    mutate(Topic = paste("Topic", as.character(Topic)), tissueSiteDetail = as.character(tissueSiteDetail), value = normalized_score) %>%
    dplyr::select(Topic, tissueSiteDetail, value)
  
  topics <- unique(sankey_data$Topic); tissues <- unique(sankey_data$tissueSiteDetail)
  topic_colors <- viridis(length(topics)); tissue_colors <- viridis(length(tissues), option = "C")
  names(topic_colors) <- topics; names(tissue_colors) <- tissues
  all_colors <- c(topic_colors, tissue_colors)
  
  ggplot(sankey_data, aes(y = value, axis1 = Topic, axis2 = tissueSiteDetail)) +
    geom_alluvium(aes(fill = Topic), alpha = alpha, width = 1/8, curve_type = "sigmoid") +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/8, color = "white") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, label.padding = unit(0.15, "lines")) +
    scale_x_discrete(limits = c("Topic", "Tissue"), expand = c(0.05, 0.05)) +
    scale_fill_manual(values = all_colors, breaks = c(names(topic_colors), names(tissue_colors))) +
    labs(title = "Topic-Tissue Association (GTEx v8 Normalized)", y = "Relative Enrichment Score (Fold Change)") +
    theme_minimal(base_size = 14) + theme(legend.position = "none", panel.grid = element_blank(), axis.text.y = element_text(size = 10))
}

create_grouped_bar_plot <- function(topic_tissue_data) {
  bar_data <- topic_tissue_data$top_tissues %>%
    mutate(Topic = paste("Topic", as.character(Topic))) %>%
    arrange(Topic, desc(normalized_score))
  
  topics <- unique(bar_data$Topic)
  topic_colors <- viridis(length(topics))
  names(topic_colors) <- topics
  
  ggplot(bar_data, aes(x = reorder(tissueSiteDetail, normalized_score), y = normalized_score, fill = Topic)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "white") +
    scale_fill_manual(values = topic_colors) +
    labs(title = "Tissue Specificity by Topic", x = "Tissue", y = "Relative Enrichment Score (Fold Change)", fill = "Topic") +
    theme_minimal(base_size = 20) + 
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
          axis.text.x = element_text(size = 18, color="black"), axis.text.y = element_text(size = 18, color="black"), panel.grid.major.y = element_blank()) +
    coord_flip()
}

create_tissue_bar_plot <- function(topic_tissue_data, top_n_tissues = 5) {
  
  plot_base_data <- topic_tissue_data$full_data %>%
    mutate(TopicNum = as.numeric(stringr::str_extract(as.character(Topic), "\\d+")))
  
  all_topics <- sort(unique(plot_base_data$TopicNum))
  
  if (length(all_topics) == 0 || all(is.na(all_topics))) {
    cat("  [Warning] No valid topic numbers extracted. Skipping faceted bar plot.\n")
    return(NULL)
  }
  
  global_palette <- scales::hue_pal()(length(all_topics))
  names(global_palette) <- as.character(all_topics)
  
  chunk_size <- 16
  total_pages <- ceiling(length(all_topics) / chunk_size)
  plot_list <- list()
  
  for (i in 1:total_pages) {
    current_topics <- all_topics[((i - 1) * chunk_size + 1):min(i * chunk_size, length(all_topics))]
    final_ncol <- min(4, length(current_topics))
    
    plot_data <- plot_base_data %>%
      filter(TopicNum %in% current_topics) %>%
      mutate(TopicNum = factor(TopicNum, levels = current_topics)) %>%
      group_by(TopicNum) %>% 
      slice_max(order_by = normalized_score, n = top_n_tissues, with_ties = FALSE) %>% 
      ungroup() %>%
      # 将重新排序的度量标准换为 normalized_score
      mutate(tissue_ordered = tidytext::reorder_within(str_wrap(tissueSiteDetail, 25), normalized_score, TopicNum))
    
    # 绘图指标全面指向 normalized_score
    p <- ggplot(plot_data, aes(x = tissue_ordered, y = normalized_score, fill = as.character(TopicNum))) +
      geom_col(alpha = 0.9, width = 0.75) + 
      scale_fill_manual(values = global_palette) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) + 
      facet_wrap(~ TopicNum, scales = "free", ncol = final_ncol) + 
      coord_flip() +
      tidytext::scale_x_reordered() + 
      labs(
        title = "Top Tissues by Topic Expression", 
        x = NULL, 
        y = "Relative Enrichment Score (Fold Change)" 
      ) +
      theme_minimal(base_size = 15) + 
      theme(
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 22, margin = margin(b = 20)),
        strip.background = element_rect(fill = "#f2f2f2", color = NA),
        strip.text = element_text(face = "bold", size = 18, color = "black", margin = margin(t=8, b=8)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(2.0, "lines"),
        panel.spacing.y = unit(1.5, "lines"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, color = "gray30"),
        axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 15)),
        axis.ticks.x = element_line(color = "gray70"), 
        axis.ticks.y = element_blank()                  
      )
    
    plot_list[[i]] <- p
  }
  return(plot_list)
}

#' Main function: Run complete tissue expression analysis (STRICTLY ONLINE REAL DATA)
run_tissue_expression_analysis <- function(top_genes_df, output_dir, top_tissues_per_topic = 5, ...) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(tidytext))
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("======================================================\n")
  cat("Starting GTEx v8 Tissue Analysis (STRICTLY ONLINE REAL DATA)\n")
  cat("======================================================\n")
  
  if (is.matrix(top_genes_df)) {
    top_genes_df <- as.data.frame(top_genes_df, stringsAsFactors = FALSE)
  }
  
  if (is.data.frame(top_genes_df)) {
    if (!"Topic" %in% colnames(top_genes_df) || !"Comorbidity" %in% colnames(top_genes_df)) {
      cat("  [Info] Converting wide top_genes_df to long format...\n")
      
      if (all(sapply(top_genes_df, is.numeric))) {
        cat("  [Info] Numeric weight matrix detected (Gaussian Style).\n")
        
        if (ncol(top_genes_df) >= nrow(top_genes_df)) {
          r_names <- rownames(top_genes_df)
          if (!any(grepl("Topic", r_names, ignore.case = TRUE))) r_names <- paste("Topic", r_names)
          top_genes_df$Topic <- r_names
          top_genes_df <- tidyr::pivot_longer(top_genes_df, cols = -Topic, names_to = "Comorbidity", values_to = "Weight")
        } else {
          c_names <- colnames(top_genes_df)
          if (!any(grepl("Topic", c_names, ignore.case = TRUE))) colnames(top_genes_df) <- paste("Topic", c_names)
          top_genes_df$Comorbidity <- rownames(top_genes_df)
          top_genes_df <- tidyr::pivot_longer(top_genes_df, cols = -Comorbidity, names_to = "Topic", values_to = "Weight")
        }
        
        cat("  [Info] Extracting Top 50 driving genes per topic for GTEx analysis...\n")
        top_genes_df <- top_genes_df %>%
          dplyr::group_by(Topic) %>%
          dplyr::slice_max(order_by = Weight, n = 50, with_ties = FALSE) %>%
          dplyr::ungroup()
        
      } 
      else if (all(sapply(top_genes_df, is.character)) || all(sapply(top_genes_df, is.factor))) {
        cat("  [Info] Character matrix detected (LDA Style).\n")
        top_genes_df <- tidyr::pivot_longer(
          top_genes_df, 
          cols = everything(), 
          names_to = "Topic", 
          values_to = "Comorbidity"
        )
      } 
      else {
        cat("  [Info] Attempting to auto-fix column names...\n")
        if (ncol(top_genes_df) >= 2) {
          colnames(top_genes_df)[1:2] <- c("Topic", "Comorbidity")
        }
      }
    }
  } else if (is.list(top_genes_df)) {
    top_genes_df <- do.call(rbind, lapply(names(top_genes_df), function(n) {
      data.frame(Topic = n, Comorbidity = top_genes_df[[n]], stringsAsFactors = FALSE)
    }))
  } else {
    stop("Error: top_genes_df format is not supported (Expected data.frame, matrix, or list).")
  }
  
  if (!is.factor(top_genes_df$Topic)) top_genes_df$Topic <- as.factor(top_genes_df$Topic)
  unique_genes <- extract_unique_genes(top_genes_df$Comorbidity)
  
  gene_mapping <- validate_and_map_genes(unique_genes)
  
  if (length(gene_mapping) == 0) {
    stop("\n[CRITICAL ERROR] Failed to map ANY genes to GTEx ENSEMBL IDs! \nPlease check your internet connection or GTEx API status. Fake data is NO LONGER permitted.")
  }
  
  expression_data <- get_gtex_expression(gene_mapping)
  
  if (is.null(expression_data) || nrow(expression_data) == 0) {
    stop("\n[CRITICAL ERROR] Downloaded expression data is empty! Halting analysis.")
  }
  
  cat("\nAnalyzing topic-tissue associations...\n")
  topic_tissue_data <- analyze_topic_tissue_association(top_genes_df, expression_data, top_tissues_per_topic)
  
  if (nrow(topic_tissue_data$full_data) == 0) return(NULL)
  
  topic_tissue_data$full_data <- prepare_visualization_data(topic_tissue_data)
  
  write.csv(topic_tissue_data$full_data, file.path(output_dir, "gtex_tissue_enrichment_full.csv"), row.names = FALSE)
  write.csv(topic_tissue_data$top_tissues, file.path(output_dir, "gtex_tissue_enrichment_top.csv"), row.names = FALSE)
  
  sankey_data <- prepare_sankey_data(topic_tissue_data)
  if (!is.null(sankey_data) && nrow(sankey_data$links) > 0) {
    saveNetwork(create_sankey_plot(sankey_data), file = file.path(output_dir, "topic_tissue_gtex_interactive.html"))
  }
  
  static_sankey <- create_static_sankey_plot(topic_tissue_data)
  if (!is.null(static_sankey)) ggsave(file.path(output_dir, "topic_tissue_gtex_sankey.png"), static_sankey, width = 16, height = 12, bg="white")
  
  bar_plot <- create_grouped_bar_plot(topic_tissue_data)
  if (!is.null(bar_plot)) ggsave(file.path(output_dir, "topic_tissue_gtex_barplot.png"), bar_plot, width = 14, height = 10, bg="white")
  
  bar_plots_list <- create_tissue_bar_plot(topic_tissue_data, top_n_tissues = top_tissues_per_topic)
  if (!is.null(bar_plots_list)) {
    for(i in seq_along(bar_plots_list)) {
      ggsave(file.path(output_dir, paste0("topic_tissue_gtex_faceted_bar_part", i, ".png")), bar_plots_list[[i]], width=16, height=12, bg="white")
    }
  }
  
  cat("\nGTEx Tissue Analysis completed successfully with REAL data!\n")
  return(list(data = topic_tissue_data, static_sankey = static_sankey, bar_plots_list = bar_plots_list))
}