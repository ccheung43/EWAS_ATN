# Creating QQ-plot for each biomarker after inflation correction

# Load required libraries 
library(tidyverse)
library(CMplot)

# Define biomarkers and subgroups
biomarkers <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
subgroups <- c("E41", "E40")

# Read in CPG annotation data 
cpg_annotation <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7)

# Function to process each biomarker-subgroup pair
process_qq <- function(biomarker, subgroup) {
  setwd("C:/Users/caitl/OneDrive - Beth Israel Lahey Health")
  
  # Read in the inflation-corrected EWAS data for each subgroup
  input_file <- paste0("2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified/", biomarker, "_", subgroup, "_EWAS_Results_Corrected.csv")
  
  df <- read.csv(input_file)
  
  group_name = paste0(biomarker, "_", subgroup, "_bc")
  biomarker_title <- switch(biomarker,
                            "ABETA4240" = paste0("A", "\u03B2", "42/40"),
                            "PTAU181" = "pTau-181",
                            "NFLIGHT" = "NfL", 
                            "GFAP" = "GFAP")
  E4_title <- switch(subgroup, 
                      "E41" = paste0("\u03B5", "4 Carrier"), 
                      "E40" = "Control")
  title = paste0(biomarker_title, " - ", E4_title)
  
  
  # Merge with CpG annotation data
  pvals <- cpg_annotation %>%
    mutate(CpG = Name) %>%
    inner_join(df, by = "CpG") %>%
    mutate(!!group_name := pval.bacon) %>%
    select(SNP = CpG, Chromosome = CHR, Position = MAPINFO, !!group_name)
  
  # FDR-adjusted p-values
  FDR_pvals <- p.adjust(pvals[[group_name]], method = "BH")
  
  # Identify significant SNPs
  SNPs_10 <- pvals$SNP[FDR_pvals < 0.10]
  SNPs_5 <- pvals$SNP[FDR_pvals < 0.05]
  
  # Determine FDR thresholds
  threshold_10 <- if (length(SNPs_10) > 0) max(pvals[[group_name]][FDR_pvals < 0.10]) else NULL
  threshold_5 <- if (length(SNPs_5) > 0) max(pvals[[group_name]][FDR_pvals < 0.05]) else NULL
  
  
  # Save Manhattan plot
  output_path <- paste0("2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified", biomarker, "_", subgroup, "_ManhattanPlot.jpg")
  
  thresh_lty <- if (is.null(threshold_5)) 2 else c(1, 2)
  
  setwd("C:/Users/caitl/OneDrive - Beth Israel Lahey Health/2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified")
  CMplot(pvals %>% select(SNP, Chromosome, Position, !!sym(group_name)), plot.type="q",col="black", multracks=FALSE, 
         box=FALSE,file="jpg",file.name="",dpi=300,
         conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
         file.output=TRUE,verbose=TRUE,width=5,height=5, main = title)
  
}

# Apply function using nested lapply()
lapply(biomarkers, function(biomarker) {
  lapply(subgroups, function(subgroup) {
    process_qq(biomarker, subgroup)
  })
})

