# SNPic: SNP Topic Modeling for Interpretable Clustering of Complex Phenotypes

We present SNPic, a probabilistic framework that redefines the analytical landscape of complex trait genetics. By conceptualizing genetic associations as a highly structured probabilistic language, SNPic deconstructs fragmented, biobank-scale GWAS catalogs into an interpretable lexicon of ``genetic topics''. These inferred topics serve as fundamental, reusable biological modules that successfully dismantle rigid clinical boundaries, exposing the true pleiotropic spectrum of human pathology. Validated through rigorous mathematical simulations, stability-optimized inference on massive human cohorts, and generalization across diverse plant and animal species, our findings demonstrate that a generative, mixed-membership approach is essential for capturing the interconnected reality of the genome. Ultimately, SNPic shifts the field from cataloging isolated variants toward reconstructing an interpretable knowledge graph of the genome's latent semantic architecture. By providing a highly scalable, privacy-preserving, and biologically transparent analytical lens, SNPic establishes a powerful new cornerstone for integrative genomics, paving the way for next-generation patient stratification and precision medicine. This work also builds a conceptual bridge between two traditionally separate fields: statistical genetics and NLP, suggesting that advances in one domain may directly transfer methodological innovations to the other.

---

## 📂 Repository Structure

The core analysis pipeline has been highly modularized into a unified dispatcher architecture:

* `run_snpic.R`: The **MAIN CODE** (CLI script) for running the SNPic pipeline. It dispatches tasks based on your chosen algorithm and matrix mode.   
* `snpic_ss_dependency.R`: Manages environment setup and dynamically loads all necessary sub-modules.   
* **Downstream Dispatchers:**   
  * `snpic_geneasword_downstream_lda.R` & `snpic_geneasword_downstream_gaussian.R`: Handles pathway, tissue, and network analyses for the **Gene-as-Word** mode.
  * `snpic_sumstat_asword_downstream_lda.R` & `snpic_sumstat_asword_downstream_gaussian.R`: Handles network and topic-trait mappings for the **Sumstat-as-Word** mode.
* **Core Models**:
  * `LDA_for_snpic_ss.R`: Implements the Latent Dirichlet Allocation (LDA) based SNPic model.
  * `mixed_membership_topics.R`: Implements the Gaussian Mixed-Membership Topic Model.
* **Downstream Analysis**:
  * `snpic_ss_similarity_analysis.R`: Calculates and visualizes disease similarities and network graphs.
  * `snpic_ss_pathway_enrichment_analysis.R`: Performs GO/KEGG pathway enrichment on genetic topics.
  * `snpic_ss_human_protein_analysis.R`: Analyzes tissue-specific expression using GTEx API.
  * `snpic_ss_ground_truth_comparison.R`: Evaluates the model against ground-truth benchmarks in simulated genotypes from LDAK.

---

## ⚙️ Dependencies & Installation

### Prerequisites
This project requires **R (>= 4.5)**. 

The analysis pipeline relies on a variety of R packages for data manipulation, topic modeling, dimensionality reduction, network visualization, and pathway enrichment. 

* **CRAN Packages**: `tidyverse`, `data.table`, `topicmodels`, `umap`, `igraph`, `ggplot2`, `optparse`, etc.
* **Bioconductor Packages**: `clusterProfiler`, `org.Hs.eg.db`

### Option 1: Conda Environment (Highly Recommended)
For local or server users, we highly recommend using Conda. This method automatically handles all complex underlying C++ libraries (e.g., `gsl` for `topicmodels`, `libxml2` for `tidyverse`) and safely avoids compilation errors.

Simply run the following commands in your terminal using the provided `environment.yml` file:

```bash
# 1. Create the environment from the provided YAML file
conda env create -f environment.yml

# 2. Activate the environment
conda activate snpic_env
```

### Option 2: Manual R Console Installation
If you prefer to install packages manually or do not have Conda installed, you can run the following script in your R console. It will automatically detect missing packages and install them:

```r
# 1. Define required CRAN packages
cran_packages <- c(
  "data.table", "dplyr", "tidyr", "reshape2", "stringr", "readr", "readxl", "tidyverse", "tidytext", "abind",
  "optparse", "httr", "jsonlite",
  "topicmodels", "tm", "maptpx", "umap", "dendextend", "proxy", "seriation", "energy", "dimRed",
  "ggplot2", "ggrepel", "ggpubr", "ggraph", "igraph", "UpSetR", "pheatmap", 
  "plotly", "networkD3", "ggalluvial", "RColorBrewer", "viridis", "scales"
)

# Install missing CRAN packages
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

# 2. Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")  
}
bioc_packages <- c("clusterProfiler", "org.Hs.eg.db")
new_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[,"Package"])]
if(length(new_bioc_packages)) BiocManager::install(new_bioc_packages)

print("All dependencies have been successfully installed!")
```
**Fallback for `clusterProfiler` Installation Errors:**   
```bash
conda install -c conda-forge -c bioconda bioconductor-clusterprofiler
```


> **📝 Notes on Installation:**
> * **`maptpx`**: Used for mixed-membership topic modeling.
> * **`clusterProfiler`**: Requires Bioconductor. Depending on your system, it might prompt you to update other Bioconductor dependencies. It is generally recommended to update them (press `a` for all) to avoid version conflicts.
> * **⚠️ Installation Errors**: If you choose Option 2 and encounter C++ compilation errors when installing `clusterProfiler` or `topicmodels`, please switch to **Option 1 (Conda)**, as it resolves these system-level dependencies seamlessly.
If you encounter C++ compilation errors when installing `clusterProfiler` natively in R, we highly recommend bypassing the scratch compilation and installing it directly using Conda. You can resolve this by running the following command in your terminal:
> * **⏳ Time Expectation**: Installation may take around **20-30 minutes** depending on your network and environment due to the large number of bioinformatics dependencies.

## 📊 Input Data Specifications 

While we provide demo data to help you get started quickly, **we highly encourage you to apply SNPic to your own custom datasets!** You can always open the files in our `data/`  and `data/` `code/snpic_ss_definitive_edition/` directory to inspect the exact formatting. 

To run your own analysis, you need to prepare three main types of input files:

### 1. GWAS Significant SNP Lists (`--input_folder`)
A directory containing a set of text files (typically `.list` or `.txt`), where each file corresponds to a specific disease or trait. These files contain the significant SNPs you have pre-filtered (e.g., $P < 5 \times 10^{-8}$).
* **Format:** Tab or space-separated.
* **Columns:**
  * **Column 1 (Required):** SNP rsID (e.g., `rs123456`).
  * **Column 2 (Optional):** P-value or summary statistics. *(Note: SNPic's core model only requires the list of variant IDs, so the second column is entirely optional and will be ignored if present).*

| rsid | pval (Optional) |
| :--- | :--- |
| rs123456 | 5e-9 |
| rs987654 | 1e-12 |

### 2. SNP-to-Gene Mapping File (`--snp_gene_map`)
This file is required if you are using the **Gene-as-Word** mode (`--mode "gene"`). It projects the genetic variants onto biological functional units (genes).
* **Format:** Tab or comma-separated.
* **Columns:** Exactly two columns.
  * **Column 1:** SNP rsID.
  * **Column 2:** Corresponding Gene Name.

| rsid | gene |
| :--- | :--- |
| rs123456 | SLC2A9 |
| rs987654 | HLA-DQA1 |

### 3. Trait Meta-Information File (`--master_map`)
A `.csv` file that links your input files to their biological annotations. 
> ⚠️ **Crucial Note:** This file **does NOT participate in the underlying LDA mathematical inference**. It is strictly used by the downstream modules to make your final network and scatter plots more intuitive, aesthetically pleasing, and scientifically grouped.

* **Format:** Comma-separated (`.csv`).
* **Columns:** Exactly four columns with the following headers:
  * `filename`: The exact name of the file in your input folder (e.g., `gout_finngen.list`).
  * `trait`: The user-friendly, display name of the trait you want to appear in the final output plots (e.g., `Gout`).
  * `label1`: A categorical variable used to assign **colors** in the downstream plots (e.g., grouping by disease system: `Metabolic`, `Autoimmune`).
  * `label2`: A categorical variable used to assign **node shapes** in the network plots (e.g., grouping by data source: `FinnGen`, `UKBB`).

**Example `master_disease_mapping.csv`:**
| filename | trait | label1 | label2 |
| :--- | :--- | :--- | :--- |
| gout_finngen.list | Gout | Metabolic | FinnGen |
| ad_ukbb.list | Alzheimer's Disease | Psychiatric | UKBB |
| sle_finngen.list | Systemic Lupus | Autoimmune | FinnGen |

## 💻 How to run SNPic

The SNPic pipeline utilizes a single, powerful global dispatcher.
First, give the script execution permissions:   
```
chmod +x run_snpic.R
```      
Then you can run the script by passing the necessary arguments using ***Rscript***.   
You can view the full list of options at any time by running:   

```
Rscript run_snpic.R --help
```
### Key Routing Parameters:
* `--mode`: Choose your base matrix construction. Options: `"gene"` (Gene-as-word) or `"ss"` (Sumstat-as-word). **Default: `"gene"`**
* `--model`: Choose the downstream probabilistic algorithm. Options: `"lda"`, `"gaussian"`, or `"both"`. **Default: `"lda"`**   
* 💡 Pro Tip: Add the `--keep_all_traits` flag (no arguments needed) to skip the stability threshold filter and keep all traits for downstream analysis and network plotting. **Default: filter out unstable traits.**

### Example A: Gene-as-word (LDA)
***Note 1*:** `--snp_gene_map` is strictly required when using `--mode "gene"`.   
***Note 2:*** You can set `--keep_all_traits` to keep all traits unfiltered.    
***Note 3:*** The author suggests you to start with **LDA** than **Gaussian**, would be faster and more informative.   
```R
Rscript run_snpic.R \
  --input_folder "<YOUR_CUSTOM_PATH>/" \
  --snp_gene_map "<YOUR_CUSTOM_PATH>/snp_gene_map_merged_finngen_ukbb.txt" \
  --master_map "<YOUR_CUSTOM_PATH>/master_disease_mapping.csv" \
  --out_prefix "<YOUR_CUSTOM_PATH>/<PREFIX>" \
  --k_min 5 \
  --k_max 20 \
  --n_bootstrap 50 \
  --cores 16 \
  --mode "gene" \
  --model "lda"
```
  
### Example B: Sumstat-as-word (Both LDA & Gaussian)
***Note 1*:** `--snp_gene_map` is **NOT** required when using `--mode "ss"`.   
***Note 2:*** You can set `--keep_all_traits` to keep all traits unfiltered.   
```R
Rscript run_snpic.R \
  --input_folder "<YOUR_CUSTOM_PATH>/sig_snp_list/" \
  --master_map "<YOUR_CUSTOM_PATH>/master_disease_mapping.csv" \
  --out_prefix "<YOUR_CUSTOM_PATH>/ss_as_word_results" \
  --k_min 5 \
  --k_max 20 \
  --n_bootstrap 50 \
  --cores 16 \
  --mode "ss" \
  --model "both" \
  --keep_all_traits
```

### 🛠️ Data Preparation: SNP-to-Gene Mapping

If you are using the **Gene-as-Word** approach, you must first prepare a SNP-to-Gene mapping list. You can choose from the following three options:

#### Option A: Use the Built-in Mapping File (Recommended)
We provide a fully processed, deduplicated, and protein-coding-filtered mapping file directly within this repository. This file contains approximately **400,000 SNPs** representing the intersection of the FinnGen and UK Biobank (UKBB) cohorts, strictly mapped to GTEx protein-coding regions. 

You can directly reference this file via the `--snp_gene_map` argument:
`./SNPic/code/snpic_ss_definitive_edition/snp_gene_map_merged_finngen_ukbb.txt`

#### Option B: Download via Google Drive
If you prefer to download the mapping file externally, you can access the exact same file using the link below:
* 📥 [Download `snp_gene_map_merged_finngen_ukbb.txt` from Google Drive](https://drive.google.com/file/d/1gWYI8Ut0_4SX-3plp-smkTlw34jt6onL/view?usp=drive_link)

Simply place this downloaded file into your `data/` or working directory and reference it via the `--snp_gene_map` argument.

#### Option C: Provide Your Own Personalized Mapping
If you wish to use your own dataset for the SNP-to-Gene projection, you can input a custom file. The file must be strictly formatted with **two columns**:
1. **Column 1:** SNP identifier
2. **Column 2:** Gene name

*(Optional utility: If you are generating this mapping from scratch, we provide an automated R script `SNPic/code/other/gene_noncoding_remove.R`. This script connects to the HGNC database to rigorously filter out non-coding RNAs and pseudogenes from your raw lists, ensuring high biological relevance).*

---

## 🚀 Quick Start with Demo Data

We provide demo data in the `data/sig_snp_list` directory so you can test the pipeline immediately after cloning the repository. The information for labelling these data is stored in `code/snpic_ss_definitive_edition/master_disease_mapping.csv`, please input it via the `--master_map` command.

1. Clone the repository and navigate to the directory:
```bash
git clone https://github.com/KirakiraZLY/SNPic.git   
cd SNPic
```
2. Make the main script executable:
```
chmod +x run_snpic.R
```
3. Run the pipeline using the provided demo data (Gene-as-word + LDA). **You can copy and paste this directly**:
```
Rscript ./run_snpic.R \
  --input_folder "./data/sig_snp_list/" \
  --snp_gene_map "./code/snpic_ss_definitive_edition/snp_gene_map_merged_finngen_ukbb.txt" \
  --master_map "./code/snpic_ss_definitive_edition/master_disease_mapping.csv" \
  --out_prefix "./results/demo_run/snpic_res" \
  --k_min 5 \
  --k_max 40 \
  --n_bootstrap 50 \
  --cores 16 \
  --mode "gene" \
  --model "lda"
```
(Note: Please ensure the mock paths `--snp_gene_map` and `--master_map` point to the correct files in your cloned `data/` folder).

### 📊 Expected Outputs
Upon successful execution, the specified out_prefix directory will contain:

> 1. Model Selection & Stability: Dual-axis K-selection plots and bootstrap stability evaluation matrices.
> 2. Topic distribution and Top-gene per topic.
> 3. Similarity Networks: Thresholded network graphs highlighting shared genetic architectures among traits.
> 4. Downstream Annotations:
> > + GO/KEGG pathway enrichment dot plots and detailed table.
> > + Topic-Tissue expression associations mapped via GTEx.


## 🎒 Our paper
doi: please wait a bit 🥹🥹🥹   
Enjoy reading. 😈😈😈
![SNPic description](./figure/fig1_graphical_abstract.png)

## ⚖️ License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 

The reproducible computational environment and associated code are also permanently archived on Code Ocean and are distributed under the same open-source license.