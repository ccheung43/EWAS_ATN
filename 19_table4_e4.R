# Creating table 4 of significant cpg sites, stratified by e4 status 

# Load required libraries 
library(tidyverse)
library(openxlsx)


# Define biomarkers and subgroups
biomarkers <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
subgroups <- c("E41", "E40")

# Read in CPG annotation data 
cpg_annotation <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7) %>% 
  mutate(CpG = Name)

# Function to process each biomarker-subgroup pair
find_sig_sites <- function(biomarker) {
  setwd("C:/Users/caitl/OneDrive - Beth Israel Lahey Health")
  
  # Read in the inflation-corrected EWAS data for each subgroup
  E41 <- read.csv(paste0("2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified/", biomarker, "_E41_EWAS_Results_Corrected.csv")) %>% mutate(Subgroup = "E41")
  E40 <- read.csv(paste0("2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified/", biomarker, "_E40_EWAS_Results_Corrected.csv")) %>% mutate(Subgroup = "E40")
  
  df <- rbind(E41, E40) %>% 
    mutate(Biomarker = biomarker) %>% 
    mutate(Biomarker = recode(Biomarker, 
                              "ABETA4240" = paste0("A", "\u03B2","42/40"), 
                              "PTAU181" = "pTau-181", 
                              "NFLIGHT" = "NfL", 
                              "GFAP" = "GFAP")) %>% 
    mutate(Subgroup = recode(Subgroup, 
                             "E41" = "E4", 
                             "E40" = "Control")) %>% 
    inner_join(cpg_annotation, by = "CpG") %>% 
    mutate(CpG_site = CpG, Chr = Chromosome_36, Nearest_gene = UCSC_RefGene_Name, 
           Effect_estimate_beta = coef, SE = sd, p_value = pval.bacon) %>% 
    select(Subgroup, Biomarker, CpG_site, Chr, Nearest_gene, Effect_estimate_beta, SE, p_value)
  
  
  # FDR-adjusted p-values
  FDR_pvals <- p.adjust(df$p_value, method = "BH")
  
  # Identify significant SNPs
  sig_sites <- df$CpG_site[FDR_pvals < 0.10]
  sig_df <- df %>% filter(CpG_site %in% sig_sites) %>% 
    arrange(p_value) %>% 
    pivot_wider(names_from = "Subgroup", values_from = c("Effect_estimate_beta", "SE", "p_value"))
  
  #top <- df %>% arrange(p_value) %>% slice_head(n = 10)
  return(sig_df)  
}

# Apply function using nested lapply()
sig_sites <- lapply(biomarkers, function(biomarker) {
  find_sig_sites(biomarker)
})

significant_sites <- do.call(rbind, sig_sites) %>%
  select("Biomarker", "CpG_site","Chr", "Nearest_gene", "Effect_estimate_beta_E4", "SE_E4", "p_value_E4", "Effect_estimate_beta_Control", "SE_Control", "p_value_Control")

row1 = c("", "", "", "", paste0("\u03B5", "4 Carrier"), "","", "Control", "", "")
row2 = c("Biomarker", "CpG site","Chr", "Nearest gene", "Effect estimate (beta)", "SE", "p-value", "Effect estimate (beta)", "SE", "p-value")

significant_sites <- rbind(row1, row2, significant_sites)

write.xlsx(significant_sites, "2024_EWAS_ATN/ATN_EWAS/Results/E4_stratified/20250507_table4_e4_significant_sites.xlsx", rowNames = FALSE, colNames = FALSE)
