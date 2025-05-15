# Creating a table for Table S5: Interaction of global genetic ancestry proportion and APOE alleles on ATN biomarkers. 
# by testing linear combination of regression coefficients for proportion genetic ancestry
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(readxl)
library(stringr) 
library(openxlsx)

# file names
dat_file = "2024_APOE_ATN_biomarkers/Data/20250124_genetic_ancestry_linear_combinations.xlsx"
table_file = "2024_APOE_ATN_biomarkers/Results/20250328_table_s5.xlsx"

# Read in data from xlsx file 
df <- read_excel(dat_file)

df <- df %>% 
  select(biomarker, allele, ancestry, proportion, est, paper_pval) %>% 
  pivot_wider(names_from = allele, values_from = c("est", "paper_pval")) %>% 
  select(biomarker, ancestry, proportion, est_E2, paper_pval_E2, est_E4, paper_pval_E4) %>% 
  mutate(ancestry = factor(ancestry, levels = c("AFRICAN", "AMERINDIAN", "EUROPEAN"), 
                           labels = c("African", "Amerindian", "European"))) %>% 
  arrange(ancestry) %>% 
  mutate(biomarker = factor(biomarker, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"), 
                            labels = c(paste0("A", "\u03B2","42/40"), "pTau-181", "NfL", "GFAP"))) %>%
  arrange(biomarker) %>% 
  mutate(biomarker = as.character(biomarker), ancestry = as.character(ancestry))
  
row1 <- c("","", "", paste0("APOE ", "\u03B5", "2"), "", paste0("APOE ", "\u03B5", "4"), "")
row2 <- c("Biomarker", "Genetic Ancestry", "Proportion", "est[CI]", "p-value", "est[CI]", "p-value")
df <- rbind(row1, row2, df)
  
write.xlsx(df, table_file, overwrite = TRUE, rowNames = FALSE, colNames = FALSE)
  
  
