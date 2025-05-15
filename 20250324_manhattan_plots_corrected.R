# Creating Manhattan plot for each biomarker after inflation correction

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


# Find significant SNPs based on FDR
biomarkers = c("ABETA4240_bc", "PTAU181_bc", "NFLIGHT_bc", "GFAP_bc")
SNPs_10 <- setNames(lapply(biomarkers, function(x) { 
  FDR_pvals <- p.adjust(pvals[[x]], method = "BH")
  # Filter significant results based on FDR threshold
  pvals$SNP[FDR_pvals < 0.10]
}), biomarkers)
SNPs_5 <- setNames(lapply(biomarkers, function(x) { 
  FDR_pvals <- p.adjust(pvals[[x]], method = "BH")
  # Filter significant results based on FDR threshold
  pvals$SNP[FDR_pvals < 0.05]
}), biomarkers)

# Find thresholds based on FDR 
thresholds_10 <- setNames(lapply(biomarkers, function(x) { 
  FDR_pvals <- p.adjust(pvals[[x]], method = "BH")
  # Filter significant results based on FDR threshold
  significant_pvals <- pvals[[x]][FDR_pvals < 0.10]
  if(length(significant_pvals) != 0) {
    return(max(significant_pvals))
  } else { 
    return(NULL)
  }
}), biomarkers)

thresholds_5 <- setNames(lapply(biomarkers, function(x) { 
  FDR_pvals <- p.adjust(pvals[[x]], method = "BH")
  # Filter significant results based on FDR threshold
  significant_pvals <- pvals[[x]][FDR_pvals < 0.05]
  if(length(significant_pvals) != 0) {
    return(max(significant_pvals))
  } else { 
    return(NULL) 
  } 
}), biomarkers)

titles <- data.frame(biomarker = c("ABETA4240_bc", "PTAU181_bc", "NFLIGHT_bc", "GFAP_bc"), 
                     title = c(paste0("A", "\u03B2","42/40 Corrected"), "pTau-181 Corrected", "NfL Corrected", "GFAP Corrected"))

path <- getwd()
setwd(paste0(path, "/2024_EWAS_ATN/ATN_EWAS/Results"))
lapply(biomarkers, function(x){ 
  thresh_lty = if (length(thresholds_5[[x]]) == 0) 2 else c(1,2)
  CMplot(pvals %>% select(SNP, Chromosome, Position, !!sym(x)), plot.type="m",LOG10=TRUE,col=c("grey30","grey60"),
         threshold=c(thresholds_5[[x]], thresholds_10[[x]]),threshold.lty=thresh_lty,threshold.col="black",  
         #signal.col=c("red","blue"),signal.cex=1, 
         highlight=c(SNPs_5[[x]], setdiff(SNPs_10[[x]], SNPs_5[[x]])), highlight.col =  c(rep("red", length(SNPs_5[[x]])), rep("blue", length(setdiff(SNPs_10[[x]], SNPs_5[[x]])))), 
         amplify=FALSE,file="jpg",file.name="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6, main = titles$title[titles$biomarker == x])
})
