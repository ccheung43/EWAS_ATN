# Inflation and bias correction of pvals using "bacon" R package
# Assuming working directory is OneDrive - Beth Israel Lahey Health 

# Load required libraries
library(tidyverse)
library(CMplot)
library(bacon)

# Define biomarkers
biomarkers <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
# Define subgroups 
subgroups <- c("E41", "E40")

# Function to apply bacon correction and save results
correct_bacon <- function(biomarker, subgroup) {
  # Construct file paths
  input_file <- paste0("2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified/", biomarker, "_", subgroup, "_EWAS_Results_Full.csv")
  output_file <- paste0("2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified/", biomarker, "_", subgroup, "_EWAS_Results_Corrected.csv")
  
  # Read data
  df <- read.csv(input_file)
  
  # Apply bacon correction
  bc <- bacon(teststatistics = NULL, effectsizes = df$coef, standarderrors = df$sd)
  
  # Add corrected values
  df_bc <- df %>%
    mutate(
      coef.bacon = bacon::es(bc),
      sd.bacon = bacon::se(bc),
      pval.bacon = pval(bc),
      fdr.bacon = p.adjust(pval(bc), method = "fdr")
    )
  
  # Save corrected results
  write.csv(df_bc, output_file, row.names = FALSE)
  
  message("Saved corrected results for ", biomarker)
}

lapply(biomarkers, function(biomarker) {
  lapply(subgroups, function(subgroup) {
    correct_bacon(biomarker, subgroup)
  })
})
