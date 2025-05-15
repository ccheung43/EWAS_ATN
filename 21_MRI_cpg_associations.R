library(tidyverse)
library(openxlsx)
library(arrow)
library(MethParquet)
source("2024_EWAS_ATN/ATN_EWAS/Code/create_methlist.R")
source("2024_EWAS_ATN/ATN_EWAS/Code/cpg_extract.R")

# read in the data
dat <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/20250508_MRI_normalized.csv")
sig_cpg <- read.xlsx("2024_EWAS_ATN/ATN_EWAS/Results/20250326_table2_significant_sites.xlsx") %>% 
  pull(CpG.site)

# create methylation list
cpg_annotation <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7)
meth_path <-"//its.caregroup.org/research/Sofer Lab/HCHS_SOL/Methylation/SOL_INCA_Batch_corrected/Combat_V02/"
mlist <- create_methlist(db_path = meth_path,cpg_col_db='CpG',subject_annot = dat,
                         subject_col_keep='all',cpgAnnot_col_keep=c(1:2,12:13,16),
                         cpg_annot = cpg_annotation, subject_id='Row_names', cpg_col_annot='Name', 
                         gene_col_name = 'UCSC_RefGene_Name')

# Association of significant CpGs with brain MRI measures:
tcv <- lm_ewas_outcome(mlist,trait='CEREBRUMV2_TCV',out_position = c('Gene','MAPINFO'),select_sites = sig_cpg,
                       covariates_string=c('AGE','BMI','SEX','CENTER','CD4T','NK','Bcell','Mono','Gran',
                                           'GENGROUP','PC1','PC2','PC3','PC4','PC5'))
tcb <- lm_ewas_outcome(mlist,trait='tcb.res',out_position = c('Gene','MAPINFO'),select_sites = sig_cpg,
                       covariates_string=c('AGE','BMI','SEX','CENTER','CD4T','NK','Bcell','Mono','Gran',
                                           'GENGROUP','PC1','PC2','PC3','PC4','PC5'))
wmh <- lm_ewas_outcome(mlist,trait='wmh.res',out_position = c('Gene','MAPINFO'),select_sites = sig_cpg,
                       covariates_string=c('AGE','BMI','SEX','CENTER','CD4T','NK','Bcell','Mono','Gran',
                                           'GENGROUP','PC1','PC2','PC3','PC4','PC5'))
hip <- lm_ewas_outcome(mlist,trait='hip.res',out_position = c('Gene','MAPINFO'),select_sites = sig_cpg,
                       covariates_string=c('AGE','BMI','SEX','CENTER','CD4T','NK','Bcell','Mono','Gran',
                                           'GENGROUP','PC1','PC2','PC3','PC4','PC5'))
tgv <- lm_ewas_outcome(mlist,trait='tgv.res',out_position = c('Gene','MAPINFO'),select_sites = sig_cpg,
                       covariates_string=c('AGE','BMI','SEX','CENTER','CD4T','NK','Bcell','Mono','Gran',
                                           'GENGROUP','PC1','PC2','PC3','PC4','PC5'))

write.csv(tcv, "2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_tcv_cpg_assosiations.csv", row.names = FALSE)
write.csv(tcb, "2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_tcb_cpg_assosiations.csv", row.names = FALSE)
write.csv(wmh, "2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_wmh_cpg_assosiations.csv", row.names = FALSE)
write.csv(hip, "2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_hip_cpg_assosiations.csv", row.names = FALSE)
write.csv(tgv, "2024_EWAS_ATN/ATN_EWAS/Results/MRI_secondary/20250508_tgv_cpg_assosiations.csv", row.names = FALSE)
