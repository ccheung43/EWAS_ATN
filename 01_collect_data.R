# Create a df for the complete dataset by collecting and organizing the data
# Assuming directory is OneDrive - Beth Israel Lahey Health

# Load required libraries 
library(haven)
library(tidyverse)

## Read in the data
# ATN biomarkers
atn<-read_sas("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Datasets/SOL-INCA/inca2_biomarkers_analysis_2306.sas7bdat")
# Phenotype data
phen <- read.csv("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Datasets/sol_cognitive_cov_20210408.csv")
# Methylation phenotype data
meth_phen <-read.csv("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Methylation/SOL_INCA_methylation/SOL INCA DNAm Database - Oct 2023/SOL.INCA.EPIC.clean.5383.sample.info.with.cellcount.extra.clocks.SOLID.only.csv")
# Genetic PC data
gpcs <- read.csv("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Datasets/subject_annotation_2017-09-05.csv")
# MPC data 
mpcs <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/20250212_SOL_MPC.csv")
# ID mapping date 
ids <- read_sas("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/ID_mapping/hchs_dbgap_id_mapping_200727.sas7bdat") 

# Data for Secondary Analysis 
# APOE genotype data 
apoe <- read_sas("//its.caregroup.org/research/Sofer Lab/HCHS_SOL/APOE/sol_inca_apoe_inv2.sas7bdat")

ids <- ids %>% mutate(SOL_ID = subject_id, HCHS_ID = ID)
atn <- atn %>% mutate(HCHS_ID = ID) %>% select(-c("CENTER"))
phen <- phen %>% mutate(HCHS_ID = as.character(ID))
meth_phen <- meth_phen %>% mutate(METH_PLATE_NUM = Sample_Plate, METH_PLATE_POS = Sample_Well)
gpcs <- gpcs %>% mutate(SOL_ID = SUBJECT_ID, HCHS_ID = as.character(HCHS_ID))
apoe <- apoe %>% mutate(APOE_GENOTYPE = APOE_E2_E3_E4)

raw_dat <- ids %>%
  left_join(atn, by = "HCHS_ID") %>% 
  left_join(phen, by = "HCHS_ID") %>%
  left_join(meth_phen, by = "SOL_ID") %>% 
  left_join(gpcs, by = c("HCHS_ID", "SOL_ID")) %>% 
  left_join(mpcs, by = "Row_names") %>% 
  left_join(apoe, by = "SOL_ID")


# Variables to keep for analysis
pc_vars = c("PC1", "PC2", "PC3", "PC4", "PC5")
biomarker_vars = c("ABETA40", "ABETA42", "PTAU181", "NFLIGHT", "GFAP")
biomarker_batch_vars = c("PLATE_NUM", "PLATE_POS")
phenotype_vars = c("SEX", "AGE", "BMI", "CENTER", "Visit") 
survey_vars = c("PSU_ID", "STRAT", "WEIGHT_NORM_OVERALL_INCA")
cell_vars = c("CD4T","NK","Bcell","Mono","Gran")
secondary_analysis_vars = c("MCI", "APOE_GENOTYPE", "GENGROUP")


dat <- raw_dat %>% 
  mutate(SEX = GENDER, # rename gender --> sex
         AGE = AGE_V2 + YRS_BTWN_V2INCA, # calculate age at SOL/INCA
         GENGROUP = gengrp6
  ) %>% 
  mutate(PC1 = if_else(!is.na(EV1), EV1, MPC1), 
         PC2 = if_else(!is.na(EV2), EV2, MPC2), 
         PC3 = if_else(!is.na(EV3), EV3, MPC3), 
         PC4 = if_else(!is.na(EV4), EV4, MPC4),
         PC5 = if_else(!is.na(EV5), EV5, MPC5)) %>% 
  mutate(ABETA40 = if_else(ABETA40_CENSORED == "NC", ABETA40, NA), # remove all censored biomarker level entries
         ABETA42 = if_else(ABETA42_CENSORED == "NC", ABETA42, NA), 
         PTAU181 = if_else(PTAU181_CENSORED == "NC", PTAU181, NA), 
         NFLIGHT = if_else(NFLIGHT_CENSORED == "NC", NFLIGHT, NA), 
         GFAP = if_else(GFAP_CENSORED == "NC", GFAP, NA)) %>% 
  mutate(across(where(is.character), as.factor), STRAT = as.factor(STRAT)) %>% # convert character vars (and STRAT) to factor
  select(all_of(c("Row_names", "SOL_ID", "HCHS_ID", 
                  survey_vars, phenotype_vars, pc_vars,
                  cell_vars, biomarker_vars, biomarker_batch_vars, secondary_analysis_vars)))
 
  
filtered_dat <- dat %>% 
  filter(!is.na(WEIGHT_NORM_OVERALL_INCA), # include only individuals in SOL/INCA 
         Visit == "V02", # include only individuals with methylation data at Visit 2
         !is.na(PC1), # include only individuals with PC data
         !if_all(all_of(biomarker_vars), is.na), # include only individuals with at least one ATN biomarker data point
         !if_any(any_of(phenotype_vars), is.na) # must have all phenotype vars
         )%>% 
  select(-c(Visit))


# 2389 individuals 
write.csv(filtered_dat,"2024_EWAS_ATN/ATN_EWAS/Data/20250224_ATN_pheno.csv",row.names = FALSE)


