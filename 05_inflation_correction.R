# Inflation and bias correction of pvals using "bacon" R package
# Assuming working directory is OneDrive - Beth Israel Lahey Health 

# Load required libraries 
library(tidyverse)
library(CMplot)
library(bacon)


# Read in the EWAS data 
ABETA4240 <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/ABETA4240_EWAS_Results_Full.csv")
PTAU181 <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/PTAU181_EWAS_Results_Full.csv")
NFLIGHT <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/NFLIGHT_EWAS_Results_Full.csv")
GFAP <- read.csv("2024_EWAS_ATN/ATN_EWAS/Results/GFAP_EWAS_Results_Full.csv")


# Use bacon to correct the inflation 
bc <- bacon(teststatistics = NULL, effectsizes = ABETA4240$coef, standarderrors = ABETA4240$sd)
ABETA4240_bc <- data.frame(
  ABETA4240,
  coef.bacon = bacon::es(bc),
  sd.bacon = bacon::se(bc),
  pval.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE
)

bc <- bacon(teststatistics = NULL, effectsizes = PTAU181$coef, standarderrors = PTAU181$sd)
PTAU181_bc <- data.frame(
  PTAU181,
  coef.bacon = bacon::es(bc),
  sd.bacon = bacon::se(bc),
  pval.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE
)


bc <- bacon(teststatistics = NULL, effectsizes = NFLIGHT$coef, standarderrors = NFLIGHT$sd)
NFLIGHT_bc <- data.frame(
  NFLIGHT,
  coef.bacon = bacon::es(bc),
  sd.bacon = bacon::se(bc),
  pval.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE
)



bc <- bacon(teststatistics = NULL, effectsizes = GFAP$coef, standarderrors = GFAP$sd)
GFAP_bc <- data.frame(
  GFAP,
  coef.bacon = bacon::es(bc),
  sd.bacon = bacon::se(bc),
  pval.bacon = pval(bc),
  fdr.bacon = p.adjust(pval(bc), method = "fdr"),
  stringsAsFactors = FALSE
)

# save results 
write.csv(ABETA4240_bc, "2024_EWAS_ATN/ATN_EWAS/Results/ABETA4240_EWAS_Results_Corrected.csv")
write.csv(PTAU181_bc, "2024_EWAS_ATN/ATN_EWAS/Results/PTAU181_EWAS_Results_Corrected.csv")
write.csv(NFLIGHT_bc, "2024_EWAS_ATN/ATN_EWAS/Results/NFLIGHT_EWAS_Results_Corrected.csv")
write.csv(GFAP_bc, "2024_EWAS_ATN/ATN_EWAS/Results/GFAP_EWAS_Results_Corrected.csv")


