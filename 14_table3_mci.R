# Creating table 3 of significant cpg sites, stratified by mci status 

# Load required libraries 
library(tidyverse)
library(openxlsx)


# Define biomarkers and subgroups
biomarkers <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
subgroups <- c("MCI1", "MCI0")

# Read in CPG annotation data 
cpg_annotation <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip=7) %>% 
  mutate(CpG = Name)

# Function to process each biomarker-subgroup pair
find_sig_sites <- function(biomarker) {
  setwd("C:/Users/caitl/OneDrive - Beth Israel Lahey Health")
  
  # Read in the inflation-corrected EWAS data for each subgroup
  MCI1 <- read.csv(paste0("2024_EWAS_ATN/ATN_EWAS/Results/MCI_stratified/", biomarker, "_MCI1_EWAS_Results_Corrected.csv")) %>% mutate(Subgroup = "MCI1")
  MCI0 <- read.csv(paste0("2024_EWAS_ATN/ATN_EWAS/Results/MCI_stratified/", biomarker, "_MCI0_EWAS_Results_Corrected.csv")) %>% mutate(Subgroup = "MCI0")
  
  df <- rbind(MCI1, MCI0) %>% 
    mutate(Biomarker = biomarker) %>% 
    mutate(Biomarker = recode(Biomarker, 
                              "ABETA4240" = paste0("A", "\u03B2","42/40"), 
                              "PTAU181" = "pTau-181", 
                              "NFLIGHT" = "NfL", 
                              "GFAP" = "GFAP")) %>% 
    mutate(Subgroup = recode(Subgroup, 
                             "MCI1" = "MCI", 
                             "MCI0" = "Control")) %>% 
    inner_join(cpg_annotation, by = "CpG") %>% 
    mutate(CpG_site = CpG, Chr = Chromosome_36, Nearest_gene = UCSC_RefGene_Name, 
           Effect_estimate_beta = coef, SE = sd, p_value = pval.bacon) %>% 
    select(Subgroup, Biomarker, CpG_site, Chr, Nearest_gene, Effect_estimate_beta, SE, p_value)
             
  
  # FDR-adjusted p-values
  FDR_pvals <- p.adjust(df$p_value, method = "BH")
  
  # Identify significant SNPs
  sig_sites <- df$CpG_site[FDR_pvals < 0.10]
  sig_df <- df %>% filter(CpG_site %in% sig_sites) %>% 
    pivot_wider(names_from = "Subgroup", values_from = c("Effect_estimate_beta", "SE", "p_value")) %>% 
    arrange(p_value_MCI)
    
  #top <- df %>% arrange(p_value) %>% slice_head(n = 10)
  return(sig_df)  
}

# Apply function using nested lapply()
sig_sites <- lapply(biomarkers, function(biomarker) {
    find_sig_sites(biomarker)
  })

significant_sites <- do.call(rbind, sig_sites) %>%
  select("Biomarker", "CpG_site","Chr", "Nearest_gene", "Effect_estimate_beta_MCI", "SE_MCI", "p_value_MCI", "Effect_estimate_beta_Control", "SE_Control", "p_value_Control") %>% 
  rowwise() %>%
  mutate(Nearest_gene = paste(unique(strsplit(Nearest_gene, ";")[[1]]), collapse = ";")) %>%
  ungroup() %>% 
  as_tibble() %>% 
  mutate(Effect_estimate_beta_MCI = round(Effect_estimate_beta_MCI, 3), 
         SE_MCI = round(SE_MCI, 3), 
         p_value_MCI = signif(p_value_MCI, 3)) %>% 
  mutate(Effect_estimate_beta_Control = round(Effect_estimate_beta_Control, 3), 
         SE_Control = round(SE_Control, 3), 
         p_value_Control = signif(p_value_Control, 3))

row1 = c("", "", "", "", "MCI", "","", "Control", "", "")
row2 = c("Biomarker", "CpG site","Chr", "Nearest gene", "Effect estimate (beta)", "SE", "p-value", "Effect estimate (beta)", "SE", "p-value")

significant_sites <- rbind(row1, row2, significant_sites)

write.xlsx(significant_sites, "2024_EWAS_ATN/ATN_EWAS/Results/MCI_stratified/20250515_table3_mci_significant_sites.xlsx", rowNames = FALSE, colNames = FALSE)
