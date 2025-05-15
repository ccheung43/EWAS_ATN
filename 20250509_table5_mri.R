# Creating table 5: Association of ATN biomarkers-associated CpG sites with brain MRI measures 

# load required packages 
library(tidyverse)
library(openxlsx)

# load the data tables of significant sites
sig_cpg <- read.xlsx("2024_EWAS_ATN/ATN_EWAS/Results/20250326_table2_significant_sites.xlsx") %>% 
  select(CpG.site, Biomarker)
tcv <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_tcv_cpg_assosiations.csv") %>% 
  select(estimate, se, fdr_bh) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
tcb <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_tcb_cpg_assosiations.csv") %>% 
  select(estimate, se, fdr_bh) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
wmh <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_wmh_cpg_assosiations.csv") %>% 
  select(estimate, se, fdr_bh) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
hip <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_hip_cpg_assosiations.csv") %>% 
  select(estimate, se, fdr_bh) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
tgv <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_tgv_cpg_assosiations.csv") %>% 
  select(estimate, se, fdr_bh) %>% 
  mutate(across(where(is.numeric), ~ signif(.x, 3)))

# combine into table 
table5 <- cbind(sig_cpg, tcv, tcb, wmh, hip, tgv)
  
row1 <- c("", "", "TCV", "", "", "TCB", "", "", "WMH", "", "", "HIPPO", "", "", "TGV", "", "")
row2 <- c("CpG Site", "Biomarker", 
          "Effect estimate", "SE", "p-value",
          "Effect estimate", "SE", "p-value",
          "Effect estimate", "SE", "p-value",
          "Effect estimate", "SE", "p-value",
          "Effect estimate", "SE", "p-value")
table5 <- rbind(row1, row2, table5)

write.xlsx(table5, "2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250509_table5_mri_associations.xlsx", rowNames = FALSE, colNames = FALSE)
