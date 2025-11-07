# ==============================================================
# Colocalization Analysis Example: eQTLs (eQTLGen Study)
# ==============================================================

# Load required libraries
library(data.table)
library(dplyr)

# ----------------------------
# Step 1: Prepare GWAS Data
# ----------------------------

# Read SMR results filtered data
dat1 <- fread('../result/smr_result/smr_gene_cis-eQTLs_eQTLGen_study.txt')

# Read original GWAS data
gwas_data <- fread('../rawdata/gwas_data.tsv.gz')

# Calculate MAF (Minor Allele Frequency)
gwas_data$MAF <- ifelse(gwas_data$effect_allele_frequency < 0.5, 
                        gwas_data$effect_allele_frequency, 
                        1 - gwas_data$effect_allele_frequency)

# Check MAF distribution
hist(gwas_data$MAF)

# Remove SNPs with MAF = 0
gwas_data <- gwas_data[gwas_data$MAF > 0, ]

# Rename columns for standardization
gwas_data <- dplyr::rename(gwas_data,
                           chromosome = chromosome,
                           base_pair_location = base_pair_location,
                           p_value = p_value,
                           variant_id = rsID,
                           MAF = MAF)

# Remove duplicate SNPs (keep lowest p-value)
gwas_data <- gwas_data[order(gwas_data$p_value), ]
gwas_data <- gwas_data[!duplicated(gwas_data$variant_id), ]

# ----------------------------
# Step 2: Prepare QTL Data
# ----------------------------

# Get unique probe IDs from SMR results
prob <- unique(dat1$probeID)
write.table(prob, file = "/path/to/eQTL_prob.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ----------------------------
# Step 3: Extract QTL Data (Bash Commands)
# ----------------------------

# The following bash commands extract QTL data for specific probes
# Note: These commands should be run in a terminal, not in R

# Create directories for QTL data extraction
system("mkdir /path/to/eQTL_Coloc_Prepare")
system("mkdir /path/to/eQTL_Coloc_Prepare/Single")

# Extract QTL data for each probe in parallel
system("sort -u /path/to/eQTL_prob.txt | parallel -j 10 \"cat /path/to/cis-eQTLs_eQTLGen.txt | grep {} > /path/to/eQTL_Coloc_Prepare/Single/{}.txt\"")

# Merge all extracted QTL data
system("cat /path/to/eQTL_Coloc_Prepare/Single/*.txt > /path/to/eQTL_Coloc_Prepare/eQTLs_coloc.txt")

# Clean up temporary files
system("rm /path/to/eQTL_Coloc_Prepare/Single/*.txt")

# ----------------------------
# Step 4: Process QTL Data in R
# ----------------------------

# Read processed QTL data
qtl_data <- fread('/path/to/eQTL_Coloc_Prepare/eQTLs_coloc.txt')

# Add column names
colnames(qtl_data) <- c("SNP", "Chr", "BP", "A1", "A2", "Freq", "Probe", 
                        "Probe_Chr", "Probe_bp", "Gene", "Orientation", 
                        "b", "SE", "p")

# Calculate MAF for QTL data
qtl_data$MAF1 <- ifelse(qtl_data$Freq < 0.5, qtl_data$Freq, 1 - qtl_data$Freq)
qtl_data <- qtl_data[order(qtl_data$p), ]

# ----------------------------
# Step 5: Perform Colocalization Analysis
# ----------------------------

# Initialize results table
output_table <- data.frame()

# Loop through each probe for colocalization analysis
for (i in unique(dat1$probeID)) {
  # Get top SNP and position information
  top_snp <- dat1$topSNP[dat1$probeID == i]
  chr <- dat1$ProbeChr[dat1$probeID == i]
  bp <- dat1$Probe_bp[dat1$probeID == i]
  
  # Subset GWAS data by chromosome
  gwas_subset <- gwas_data[gwas_data$chromosome == chr, ]
  
  # Subset QTL data by chromosome and position window
  qtl_subset <- qtl_data[qtl_data$Chr == chr & 
                           qtl_data$BP > (bp - 1000000) & 
                           qtl_data$BP < (bp + 1000000) & 
                           qtl_data$Probe == i, ]
  
  # Merge GWAS and QTL data
  merged_data <- merge(gwas_subset, qtl_subset,
                       by.x = 'variant_id', by.y = 'SNP',
                       suffixes = c("_gwas", "_qtl"))
  
  # Handle p-values of 0
  merged_data$p_value[merged_data$p_value == 0] <- min(merged_data$p_value[merged_data$p_value != 0])
  merged_data$p[merged_data$p == 0] <- min(merged_data$p[merged_data$p != 0])
  
  # ----------------------------
  # Colocalization Analysis
  # ----------------------------
  library(coloc)
  
  # Set sample sizes (REPLACE WITH ACTUAL VALUES)
  ngwas <- "TOTAL_GWAS_SAMPLE_SIZE"    # Total GWAS sample size
  case <- "NUMBER_OF_CASES"            # Number of cases
  nqtl <- "QTL_SAMPLE_SIZE"            # QTL sample size
  
  # Perform colocalization analysis
  result <- coloc.abf(
    dataset1 = list(pvalues = merged_data$p_value,  # GWAS p-values
                    snp = merged_data$variant_id,
                    type = "cc", 
                    N = ngwas, 
                    s = case/ngwas,
                    MAF = merged_data$MAF),
    dataset2 = list(pvalues = merged_data$p,       # QTL p-values
                    snp = merged_data$variant_id,
                    type = "quant", 
                    N = nqtl,
                    MAF = merged_data$MAF1),
    p12 = 5e-05)  # Prior probability of colocalization
  
  # Store results
  tmp <- data.frame(ProbeID = i, 
                    PPH4 = result$summary[6],       # Posterior probability of colocalization
                    SNP.PP.H4 = max(result$results$SNP.PP.H4))
  output_table <- rbind(output_table, tmp)
  
  # Create results directory if needed
  if (!dir.exists('../result/coloc_result')) {
    dir.create('../result/coloc_result', recursive = TRUE)
  }
  
  # Save detailed results if PPH4 > 0.5
  if (result$summary[6] > 0.5) {
    # Save full colocalization results
    write.table(result$results,
                paste0("../result/coloc_result/coloc_cis-eQTLs_eQTLGen_study-", top_snp, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Create data directory if needed
    if (!dir.exists('../data/locusdata')) {
      dir.create('../data/locusdata', recursive = TRUE)
    }
    
    # Prepare data for locus comparison plots
    gwas_plot <- data.frame(rsid = merged_data$variant_id, 
                            pval = merged_data$p_value)
    qtl_plot <- data.frame(rsid = merged_data$variant_id, 
                           pval = merged_data$p)
    
    # Save plot data
    write.table(gwas_plot,
                paste0("../data/locusdata/gwas_cis-eQTLs_eQTLGen_study-", i, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(qtl_plot,
                paste0("../data/locusdata/qtl_cis-eQTLs_eQTLGen_study-", i, ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Create locus comparison plot
    # Note: locuscompare function should be defined elsewhere
    locuscompare(paste0("../data/locusdata/gwas_cis-eQTLs_eQTLGen_study-", i, ".tsv"),
                 paste0("../data/locusdata/qtl_cis-eQTLs_eQTLGen_study-", i, ".tsv"),
                 title1 = 'GWAS', 
                 title2 = 'eQTL',
                 legend = TRUE)
    
    # Save plot
    if (!dir.exists('../figure')) dir.create('../figure', recursive = TRUE)
    ggsave(paste0("../figure/locus_cis-eQTLs_eQTLGen_study-", i, ".jpeg"),
           width = 10, height = 6)
  }
}