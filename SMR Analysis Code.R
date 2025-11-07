# ==============================================================
# SMR Analysis Pipeline Script

# ==============================================================

# ----------------------------
# Section 1: Environment Setup and Configuration
# ----------------------------

#' Initialize Analysis Environment
#' 
#' Load required R packages and set analysis environment
initialize_environment <- function() {
  # Install necessary packages if not already installed
  required_packages <- c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", 
                         "forestplot", "clusterProfiler", "org.Hs.eg.db")
  new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
  
  # Load packages
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(forestplot)
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  
  # Set paths (User must update these according to their system)
  set_paths <- function() {
    list(
      smr_bin = "/path/to/smr-1.3.1-linux-x86_64/smr",
      bfile = "/path/to/g1000_eur/g1000_eur",
      mqtl_data = "/path/to/mqtl/LBC_BSGS_meta/meta/LBC_BSGS_meta_merge",
      eqtl_data = "/path/to/eqtl/cis-eQTLs_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense",
      pqtl_data = "/path/to/pqtl/Decode_2550/pqtl_p0.05_download_sucess_2550",
      probe_info = "/path/to/mqtl/LBC_BSGS_meta/meta/GPL16304_LBC_BSGS.txt",
      glist = "/path/to/glist-hg19.txt"
    )
  }
  
  paths <<- set_paths()
  
  message("Analysis environment initialized successfully")
}

# ----------------------------
# Section 2: GWAS Data Processing
# ----------------------------

#' Process GWAS Data
#' 
#' Convert raw GWAS data to SMR-compatible format
#' 
#' @param gwas_file Path to raw GWAS data file
#' @param output_file Path for processed output file
process_gwas_data <- function(gwas_file, output_file) {
  # Read data
  gwas_data <- fread(gwas_file)
  
  # Process columns and format
  processed <- gwas_data %>%
    select(
      SNP = rsids,
      A1 = alt,
      A2 = ref,
      freq = af_alt,
      b = beta,
      se = sebeta,
      p = pval
    ) %>%
    mutate(n = NA)  # Set sample size to NA if unknown
  
  # Save processed data
  fwrite(processed, file = output_file, sep = "\t", quote = FALSE)
  
  message("GWAS data processed and saved to: ", output_file)
  return(output_file)
}

# ----------------------------
# Section 3: SMR Analysis Execution
# ----------------------------

#' Run SMR Analysis
#' 
#' Execute SMR analysis command
#' 
#' @param gwas_summary Path to processed GWAS data file
#' @param beqtl_summary Path to QTL data
#' @param output_prefix Output file prefix
#' @param qtl_type Type of QTL ("mqtl", "eqtl", "pqtl")
run_smr_analysis <- function(gwas_summary, beqtl_summary, output_prefix, qtl_type = "mqtl") {
  # Build command
  cmd <- sprintf(
    "%s --bfile %s --gwas-summary %s --beqtl-summary %s --out %s --thread-num 10 --diff-freq 0.2 --diff-freq-prop 0.05 --smr-multi",
    paths$smr_bin,
    paths$bfile,
    gwas_summary,
    beqtl_summary,
    output_prefix
  )
  
  # Execute command
  system(cmd)
  
  # Verify output file
  output_file <- paste0(output_prefix, ".msmr")
  if (!file.exists(output_file)) {
    warning("SMR analysis may have failed. Output file not found: ", output_file)
    return(NULL)
  }
  
  message(qtl_type, " SMR analysis completed. Results saved to: ", output_file)
  return(output_file)
}

# ----------------------------
# Section 4: SMR Results Processing
# ----------------------------

#' Process SMR Results
#' 
#' Apply FDR correction and filter genes
#' 
#' @param smr_file Path to SMR results file
#' @param qtl_type Type of QTL ("mqtl", "eqtl", "pqtl")
#' @param gene_set Gene set to filter results
process_smr_results <- function(smr_file, qtl_type = "mqtl", gene_set = NULL) {
  # Read results
  results <- fread(smr_file)
  
  # Apply FDR correction
  results$FDR <- p.adjust(results$p_SMR, method = "BH")
  
  # Filter genes if gene set provided
  if (!is.null(gene_set)) {
    if (qtl_type == "mqtl") {
      # Add gene symbols for mQTL
      probe_info <- fread(paths$probe_info)
      results <- results %>%
        left_join(select(probe_info, ID, SYMBOL = Closest_TSS_gene_name), 
                  by = c("probeID" = "ID"))
      
      filtered <- results %>% filter(SYMBOL %in% gene_set)
    } else if (qtl_type == "eqtl") {
      # Convert gene symbols to ENSEMBL IDs for eQTL
      gene_map <- bitr(gene_set, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db)
      filtered <- results %>% filter(Gene %in% gene_map$ENSEMBL)
    } else { # pqtl
      filtered <- results %>% filter(probeID %in% gene_set)
    }
  } else {
    filtered <- results
  }
  
  # Save processed results
  output_file <- gsub("\\.msmr$", "_processed.txt", smr_file)
  fwrite(filtered, output_file, sep = "\t")
  
  message("SMR results processed and saved to: ", output_file)
  return(filtered)
}

# ----------------------------
# Section 5: Downstream SMR Analysis
# ----------------------------

#' Run Downstream SMR Analysis
#' 
#' Execute mQTL-eQTL or eQTL-pQTL analysis
#' 
#' @param exposure_probes File with exposure probe list
#' @param outcome_probes File with outcome probe list
#' @param output_prefix Output file prefix
#' @param analysis_type Analysis type ("m_e" or "e_p")
run_downstream_smr <- function(exposure_probes, outcome_probes, output_prefix, analysis_type = "m_e") {
  # Set QTL data paths based on analysis type
  if (analysis_type == "m_e") {
    beqtl1 <- paths$mqtl_data
    beqtl2 <- paths$eqtl_data
  } else { # e_p
    beqtl1 <- paths$eqtl_data
    beqtl2 <- paths$pqtl_data
  }
  
  # Build command
  cmd <- sprintf(
    "%s --bfile %s --beqtl-summary %s --beqtl-summary %s --extract-exposure-probe %s --extract-outcome-probe %s --out %s --smr-multi",
    paths$smr_bin,
    paths$bfile,
    beqtl1,
    beqtl2,
    exposure_probes,
    outcome_probes,
    output_prefix
  )
  
  # Execute command
  system(cmd)
  
  # Verify output file
  output_file <- paste0(output_prefix, ".msmr")
  if (!file.exists(output_file)) {
    warning("Downstream SMR analysis may have failed. Output file not found: ", output_file)
    return(NULL)
  }
  
  message(analysis_type, " SMR analysis completed. Results saved to: ", output_file)
  return(output_file)
}

# ----------------------------
# Section 6: Visualization Functions
# ----------------------------

#' Create Forest Plot
#' 
#' Generate forest plot from SMR results
#' 
#' @param data Processed SMR results data frame
create_forest_plot <- function(data) {
  # Prepare label data
  label <- cbind(
    c("probeID", data$probeID),
    c("Gene", data$SYMBOL),
    c("OR (95%CI)", sprintf("%.2f (%.2f-%.2f)", data$OR_SMR, data$CI_SMR_lower, data$CI_SMR_upper)),
    c("Pvalue", sprintf("%.3g", data$p_SMR)),
    c("FDR", sprintf("%.3g", data$FDR))
  )
  
  # Create forest plot
  forestplot(
    labeltext = label,
    mean = c(NA, data$OR_SMR),
    lower = c(NA, data$CI_SMR_lower),
    upper = c(NA, data$CI_SMR_upper),
    zero = 1,
    boxsize = 0.2,
    lwd.zero = 2,
    lwd.ci = 2,
    col = fpColors(box = 'black', lines = 'black', zero = 'grey'),
    hrzl_lines = list("1" = gpar(lty = 1, lwd = 1, col = "black")),
    graph.pos = 3,
    clip = c(0.5, 2),
    xticks = c(0.5, 1, 1.5, 2),
    txt_gp = fpTxtGp(
      label = gpar(fontfamily = "serif", cex = 0.8),
      ticks = gpar(fontfamily = "serif", cex = 0.8),
      xlab = gpar(fontfamily = "serif", cex = 0.8)
    )
  ) |> 
    fp_set_zebra_style("#EFEFEF")
}

#' Create Manhattan Plot
#' 

# Load custom Manhattan plot refinement function
source("/path/to/custom_functions/refine_manhattan_plot.R") 

# -------------------------------------------------------------------
# Identify Highlight Points for Visualization
# -------------------------------------------------------------------
#
# Strategy for selecting points to highlight:
# 1. Optimal: Genes showing significance in all three QTL types (mQTL, eQTL, pQTL)
# 2. Alternative: Genes passing cross-filtering criteria:
#    - Significant in both mQTL-GWAS and eQTL-GWAS analyses
#    - Significant in mQTL-eQTL analysis
#
# These points will be highlighted in the Manhattan plot

# Read SMR results for different QTL types
eqtl_results <- fread('smr_eQTLs_eQTLGen_study.txt')
mqtl_results <- fread('smr_mQTLs_study.txt')
pqtl_results <- fread('smr_pQTLs_study.txt')  # Replace with actual pQTL file name

# Find genes significant in all three QTL types
common_genes <- intersect(
  intersect(eqtl_results$SYMBOL, mqtl_results$SYMBOL),
  pqtl_results$probeID
) 

# Get probe IDs for these genes from mQTL data
highlight_probes <- mqtl_results$probeID[mqtl_results$SYMBOL %in% common_genes]

# -------------------------------------------------------------------
# Prepare Input Data for Plotting
# -------------------------------------------------------------------
#
# Required columns:
# - SYMBOL: Gene symbol (for labeling)
# - probeID: Probe identifier
# - ProbeChr: Chromosome location of probe
# - Probe_bp: Base pair position of probe
# - p_SMR_multi: P-value for SMR analysis (or FDR column if using that)

# Read SMR results (replace 'dat' with your actual results dataframe)
gwas_data <- fread('smr_results.msmr')[, .(SYMBOL, probeID, ProbeChr, Probe_bp, p_SMR_multi)]

# -------------------------------------------------------------------
# Generate Manhattan Plot
# -------------------------------------------------------------------
#
# Using custom refinement function with recommended parameters
Manhattan_refine(
  Background_data = gwas_data,         # Input data frame
  highlight_names = highlight_probes,    # Probes to highlight
  thresholds_line_color = "#E58601",    # Threshold line color (orange)
  points_color = "#A5C2A3",             # Background points color (green)
  highlight_color = "#972D15",           # Highlight points color (burgundy)
  pointsize = 1.5,                      # Size of background points
  highlight_size = 2,                   # Size of highlighted points
  label_size = 3.5,                     # Size of data labels
  force = 0.2,                          # Repulsion force between labels
  pval_names = "p_SMR_multi",           # P-value column to use
  thresholds = 0.05,                    # Significance threshold
  legend_position = c(0.9, 0.8),        # Legend position (x,y)
  legend_labels = "SMR: cis-mQTL Plasma" # Legend text (customize per QTL type)
)



#' Generate SMR Locus Plot
#' 
#' @param data Input data object (SMR results)
#' @param smr_thresh Significance threshold for SMR p-values
#' @param heidi_thresh Significance threshold for HEIDI test
#' @param plotWindow Genomic window size around locus (in kb)
#' @param max_anno_probe Maximum number of probes to annotate
#' @param smr_thresh_type P-value type to use for thresholding
#' @param smrindx Specific probes to highlight
SMRLocusPlot(
  data = SMRData,             # Input SMR results data
  smr_thresh = 0.05,          # Significance threshold for SMR p-values
  heidi_thresh = 0.05,        # Significance threshold for HEIDI test
  plotWindow = 500,           # Genomic window size around locus (in kb)
  max_anno_probe = 15,        # Maximum number of probes to annotate
  smr_thresh_type = "p_SMR_multi",  # P-value type for thresholding
  smrindx = probeid.list       # Specific probes to highlight
)

# ==============================================================
# SMR Effect Plot Configuration
# ==============================================================

# ---------------------------------------------------------------
# Generate SMR Effect Plot
# ---------------------------------------------------------------

#' Customized SMR Effect Visualization
#'
#' @param data Input data object (SMR results)
#' @param trait_name Main title for the plot
#' @param title1 X-axis title (exposure type)
#' @param title2 Y-axis title (outcome type)
#' @param Pval_title P-value type to display
SMREffectPlot(
  data = SMRData,             # Input SMR results data
  trait_name = "Effect Plot for Target_Gene",  # Main plot title
  title1 = "mQTLs",            # X-axis title (modify per data type)
  title2 = "GWAS",             # Y-axis title (typically unchanged)
  Pval_title = "p_SMR_multi"   # P-value type for display
)

# ----------------------------
# Section 7: Main Analysis Workflow
# ----------------------------

#' Execute Full SMR Analysis Workflow
#' 
#' @param raw_gwas_file Path to raw GWAS data file
#' @param gene_set Gene set of interest (optional)
main_analysis_workflow <- function(raw_gwas_file, gene_set = NULL) {
  # Initialize environment
  initialize_environment()
  
  # Step 1: Process GWAS data
  processed_gwas <- process_gwas_data(raw_gwas_file, "processed_gwas.txt")
  
  # Step 2: Run core SMR analyses
  mqtl_result <- run_smr_analysis(processed_gwas, paths$mqtl_data, "mqtl_results", "mqtl")
  eqtl_result <- run_smr_analysis(processed_gwas, paths$eqtl_data, "eqtl_results", "eqtl")
  pqtl_result <- run_smr_analysis(processed_gwas, paths$pqtl_data, "pqtl_results", "pqtl")
  
  # Step 3: Process SMR results
  mqtl_processed <- process_smr_results(mqtl_result, "mqtl", gene_set)
  eqtl_processed <- process_smr_results(eqtl_result, "eqtl", gene_set)
  pqtl_processed <- process_smr_results(pqtl_result, "pqtl", gene_set)
  
  # Step 4: Identify significant genes
  significant_genes <- identify_significant_genes(mqtl_processed, eqtl_processed, pqtl_processed)
  
  # Step 5: Run downstream analyses
  if (!is.null(significant_genes$m_e)) {
    # Create probe lists for m-e analysis
    create_probe_file(significant_genes$m_e$exposure_probes, "exposure_probes.txt")
    create_probe_file(significant_genes$m_e$outcome_probes, "outcome_probes.txt")
    
    # Run m-e analysis
    run_downstream_smr("exposure_probes.txt", "outcome_probes.txt", "m_e_results", "m_e")
  }
  
  # Step 6: Create visualizations
  # Forest plot
  create_forest_plot(mqtl_processed)
  
  # Manhattan plot
  highlight_probes <- if (!is.null(significant_genes)) significant_genes$all$probeID else NULL
  create_manhattan_plot(mqtl_processed, highlight_probes, "mqtl")
  
  message("Full SMR analysis workflow completed successfully")
}

# ----------------------------
# Section 8: Helper Functions
# ----------------------------

#' Identify Significant Genes
#' 
#' Find genes significant across all QTL types
#' 
#' @param mqtl_data mQTL processed results
#' @param eqtl_data eQTL processed results
#' @param pqtl_data pQTL processed results
identify_significant_genes <- function(mqtl_data, eqtl_data, pqtl_data) {
  # Filter significant genes (p_SMR_multi < 0.05)
  sig_mqtl <- mqtl_data %>% filter(p_SMR_multi < 0.05)
  sig_eqtl <- eqtl_data %>% filter(p_SMR_multi < 0.05)
  sig_pqtl <- pqtl_data %>% filter(p_SMR_multi < 0.05)
  
  # Find genes significant in all three analyses
  common_genes <- intersect(
    intersect(sig_mqtl$SYMBOL, sig_eqtl$SYMBOL),
    sig_pqtl$SYMBOL
  )
  
  # Prepare data for downstream analysis
  m_e_analysis <- list(
    exposure_probes = sig_mqtl %>% 
      filter(SYMBOL %in% common_genes) %>% 
      pull(probeID),
    outcome_probes = sig_eqtl %>% 
      filter(SYMBOL %in% common_genes) %>% 
      pull(probeID)
  )
  
  return(list(
    all = common_genes,
    m_e = m_e_analysis
  ))
}

#' Create Probe List File
#' 
#' @param probes Vector of probe IDs
#' @param filename Output file name
create_probe_file <- function(probes, filename) {
  writeLines(probes, filename)
  message("Probe list saved to: ", filename)
}

