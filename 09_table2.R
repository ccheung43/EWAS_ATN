# Creating table 2 of significant cpg sites

# Load required libraries 
library(tidyverse)
library(openxlsx)

# Read in the inflation-corrected EWAS data 
ABETA4240 <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/ABETA4240_EWAS_Results_Corrected.csv")
PTAU181 <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/PTAU181_EWAS_Results_Corrected.csv")
NFLIGHT <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/NFLIGHT_EWAS_Results_Corrected.csv")
GFAP <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/GFAP_EWAS_Results_Corrected.csv")

# Read in CPG annotation data 
cpg_annotation <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7) %>% 
  mutate(CpG = Name)


# Find significant cpgs based on FDR
biomarkers = list(ABETA4240 = ABETA4240, 
                  PTAU181 = PTAU181, 
                  NFLIGHT = NFLIGHT, 
                  GFAP = GFAP) 
cpgs_sig_10 <- setNames(lapply(biomarkers, function(biomarker) { 
  FDR_pvals <- p.adjust(biomarker$pval.bacon, method = "BH")
  # Filter significant results based on FDR threshold
  biomarker$CpG[FDR_pvals < 0.10]
}), names(biomarkers)) 


significant_sites <- data.frame(
  Biomarker = character(),
  CpG_site = character(),
  Chr = character(),
  Nearest_gene = character(),
  Effect_estimate_beta = numeric(),
  SE = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE  # Ensure character columns stay as characters
)

for(biomarker in names(biomarkers)) { 
  sites <- biomarkers[[biomarker]] %>% filter(CpG %in% cpgs_sig_10[[biomarker]]) %>% 
    inner_join(cpg_annotation, by = "CpG") %>% 
    mutate(Biomarker = biomarker, CpG_site = CpG, Chr = Chromosome_36, Nearest_gene = UCSC_RefGene_Name, 
           Effect_estimate_beta = coef, SE = sd, p_value = pval.bacon) %>% 
    select(Biomarker, CpG_site, Chr, Nearest_gene, Effect_estimate_beta, SE, p_value)
  significant_sites <- rbind(significant_sites, sites)
} 


top10_sites <- data.frame(
  Biomarker = character(),
  CpG_site = character(),
  Chr = character(),
  Nearest_gene = character(),
  Effect_estimate_beta = numeric(),
  SE = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE  # Ensure character columns stay as characters
)

for(biomarker in names(biomarkers)) { 
  sites <- biomarkers[[biomarker]] %>% arrange(pval.bacon) %>% slice_head(n = 10) %>% 
    inner_join(cpg_annotation, by = "CpG") %>% 
    mutate(Biomarker = biomarker, CpG_site = CpG, Chr = Chromosome_36, Nearest_gene = UCSC_RefGene_Name, 
           Effect_estimate_beta = coef, SE = sd, p_value = pval.bacon) %>% 
    select(Biomarker, CpG_site, Chr, Nearest_gene, Effect_estimate_beta, SE, p_value)
  top10_sites <- rbind(top10_sites, sites)
} 


significant_sites <- significant_sites %>% 
  mutate(Biomarker = recode(Biomarker, 
    "ABETA4240" = paste0("A", "\u03B2","42/40"), 
    "PTAU181" = "pTau-181", 
    "NFLIGHT" = "NfL", 
    "GFAP" = "GFAP")
  ) %>% 
  rowwise() %>%
  mutate(Nearest_gene = paste(unique(strsplit(Nearest_gene, ";")[[1]]), collapse = ";")) %>%
  ungroup() %>% 
  as_tibble() %>% 
  mutate(Effect_estimate_beta = round(Effect_estimate_beta, 3), 
         SE = round(SE, 3), 
         p_value = signif(p_value, 3))

top10_sites <- top10_sites %>% 
  mutate(Biomarker = recode(Biomarker, 
                            "ABETA4240" = paste0("A", "\u03B2","42/40"), 
                            "PTAU181" = "pTau-181", 
                            "NFLIGHT" = "NfL", 
                            "GFAP" = "GFAP")
  )


row1 = c("Biomarker", "CpG site","Chr", "Nearest gene", "Effect estimate (beta)", "SE", "p-value")

significant_sites <- rbind(row1, significant_sites)
top10_sites <- rbind(row1, top10_sites)

write.xlsx(significant_sites, "2024_EWAS_ATN/ATN_EWAS/Results/20250515_table2_significant_sites.xlsx", rowNames = FALSE, colNames = FALSE)
#write.xlsx(top10_sites, "2024_EWAS_ATN/ATN_EWAS/Results/20250327_table2_top10_sites.xlsx", rowNames = FALSE, colNames = FALSE)

