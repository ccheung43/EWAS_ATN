# Creating QQ-plot for each biomarker after inflation correction

# Load required libraries 
library(tidyverse)
library(CMplot)

# Read in the inflation-corrected EWAS data 
ABETA4240 <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/ABETA4240_EWAS_Results_Corrected.csv")
PTAU181 <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/PTAU181_EWAS_Results_Corrected.csv")
NFLIGHT <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/NFLIGHT_EWAS_Results_Corrected.csv")
GFAP <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/GFAP_EWAS_Results_Corrected.csv")

# Read in CPG annotation data 
cpg_annotation <- read.csv("2024_EWAS_ATN/ATN_EWAS/Data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",skip=7)

# Create df with EWAS pvals and cpg annotation data 
pvals <- cpg_annotation %>% mutate(CpG = Name) %>% 
  inner_join(ABETA4240_bc, by = "CpG") %>% mutate(ABETA4240_bc = pval.bacon) %>% select(-c(pval.bacon)) %>% 
  inner_join(PTAU181_bc, by = "CpG") %>% mutate(PTAU181_bc = pval.bacon) %>% select(-c(pval.bacon)) %>% 
  inner_join(NFLIGHT_bc, by = "CpG") %>% mutate(NFLIGHT_bc = pval.bacon) %>% select(-c(pval.bacon)) %>% 
  inner_join(GFAP_bc, by = "CpG") %>% mutate(GFAP_bc = pval.bacon) %>% select(-c(pval.bacon)) %>% 
  mutate(SNP = CpG, Chromosome = CHR, Position = MAPINFO) %>% 
  select(SNP, Chromosome, Position, ABETA4240_bc, PTAU181_bc, NFLIGHT_bc, GFAP_bc)

#colnames(pvals) = c("SNP", "Chromosome", "Position", paste0("A", "\u03B2","42/40"), "pTau-181", "NfL", "GFAP")

path <- getwd()
#setwd(paste0(path, "/2024_EWAS_ATN/ATN_EWAS/Results"))
CMplot(pvals,plot.type="q",col="black",multracks=TRUE,
       box=FALSE,file="jpg",file.name="",dpi=300,
       conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
       file.output=TRUE,verbose=TRUE,width=5,height=5)


biomarkers = c("ABETA4240_bc", "PTAU181_bc", "NFLIGHT_bc", "GFAP_bc")
titles <- data.frame(biomarker = biomarkers, 
                     title = c(paste0("A", "\u03B2","42/40"), "pTau-181", "NfL", "GFAP"))

lapply(biomarkers, function(x){ 
  CMplot(pvals %>% select(SNP, Chromosome, Position, !!sym(x)), plot.type="q",col="black", multracks=FALSE, 
         box=FALSE,file="jpg",file.name="",dpi=300,
         conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
         file.output=TRUE,verbose=TRUE,width=5,height=5, main = titles$title[titles$biomarker == x])
})
