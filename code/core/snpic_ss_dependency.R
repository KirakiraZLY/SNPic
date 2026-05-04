############################################
# 1. Load dependencies & Setup
############################################
cat("Loading dependencies...\n")

# Suppress startup messages to keep CLI output clean
suppressPackageStartupMessages({
  # Data manipulation & Base
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(stringr)
  library(readr)
  library(readxl)
  library(tidyverse)
  library(parallel)
  library(abind)
  
  # Plotting & Visualization
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(ggraph)
  library(igraph)
  library(UpSetR)
  library(pheatmap)
  library(plotly)
  library(networkD3)
  library(ggalluvial)
  library(RColorBrewer)
  library(viridis)
  library(scales)
  
  # Clustering, Topic Modeling & Dim Reduction
  library(topicmodels)
  library(tm)
  library(maptpx)
  library(umap)
  library(cluster)
  library(dendextend)
  library(proxy)
  library(seriation)
  library(energy)
  library(dimRed)
  
  # Bioinformatics
  library(clusterProfiler)
  
  # Web / API
  library(httr)
  library(jsonlite)
  library(tidytext)
})

set.seed(123)

############################################
# 2. Source Custom Scripts & Functions
############################################
cat("Sourcing custom analytical scripts dynamically...\n")

# Determine the base directory dynamically based on CLI argument
if (exists("script_dir")) {
  BASE_DIR <- script_dir
} else {
  BASE_DIR <- getwd() # Only as Fallback when running this script standalone (not via run_snpic.R)
}
PARENT_DIR <- dirname(BASE_DIR)

# Core SNPic Models (Same directory)
source(file.path(BASE_DIR, "LDA_for_snpic_ss.R"))
source(file.path(BASE_DIR, "mixed_membership_topics.R"))

# Downstream Analysis (Same directory)
source(file.path(BASE_DIR, "snpic_ss_similarity_analysis.R"))
source(file.path(BASE_DIR, "snpic_ss_ground_truth_comparison.R"))
source(file.path(BASE_DIR, "snpic_ss_pathway_enrichment_analysis.R"))
source(file.path(BASE_DIR, "snpic_ss_human_protein_analysis.R"))
# source(file.path(BASE_DIR, "label_disease_map.R"))
source(file.path(BASE_DIR, "snpic_geneasword_downstream_lda.R"))
source(file.path(BASE_DIR, "snpic_geneasword_downstream_gaussian.R"))
source(file.path(BASE_DIR, "snpic_sumstat_asword_downstream_lda.R"))
source(file.path(BASE_DIR, "snpic_sumstat_asword_downstream_gaussian.R"))


# Evaluation & Omics integration (Sibling directories in PARENT_DIR)
#source(file.path(PARENT_DIR, "snpic_evaluation", "snpic_dtm_evaluation.R"))
#source(file.path(PARENT_DIR, "snpic_x_multi_omics_gaussian", "snpic_x_read_individual_gwas_dtm.R"))
#source(file.path(PARENT_DIR, "snpic_evaluation", "snpic_interpretability_evaluation_shannon_entropy.R"))

cat("Setup complete. Starting analysis...\n")
############################################