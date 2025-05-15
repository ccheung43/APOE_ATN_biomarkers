# Organizing the data and selecting desired variables
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(dplyr)
library(tidyverse)
library(naniar)


# Read in data from csv file 
data_file = "2024_APOE_ATN_biomarkers/Data/20241001_data.csv"
full_dat <- read.csv(data_file)
plot_file = "2024_APOE_ATN_biomarkers/Results/20241028_missing_upset_plot.pdf"
missing_data_file = "2024_APOE_ATN_biomarkers/Data/20241024_data_missing.csv"
complete_data_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_complete.csv"

# Create a df to contain desired exposures, outcomes, and covariates 
exposure_vars = c("APOE_E2_E3_E4")
outcome_vars = c("ABETA40", "ABETA42", "PTAU181", "NFLIGHT", "GFAP")
outcome_censoring = c("ABETA40_CENSORED", "ABETA42_CENSORED", "PTAU181_CENSORED", "NFLIGHT_CENSORED", "GFAP_CENSORED")
basic_covars = c("GENDER", "AGE_V2", "YRS_BTWN_V2INCA", "CENTER", "EV1", "EV2", "EV3", "EV4", "EV5")
stratification_vars = c("gengrp6", "mean_afr", "mean_eur", "mean_amer")
survey_vars = c("PSU_ID", "STRAT", "WEIGHT_NORM_OVERALL_INCA")
batch_vars = c("PLATE_NUM", "PLATE_POS")

reduced_dat <- full_dat %>% 
  mutate(SOL_ID = coalesce(SOL_ID, subject_id, subject_id.x, subject_id.y, SUBJECT_ID.x, SUBJECT_ID.y)) %>% # fixing the SOL_ID column because it somehow got added into multiple different cols 
  select(all_of(c("SOL_ID", exposure_vars, outcome_vars, outcome_censoring, basic_covars, stratification_vars, survey_vars, batch_vars))) %>% 
  filter(!is.na(WEIGHT_NORM_OVERALL_INCA)) %>% 
  mutate(SEX = GENDER, 
         AGE = AGE_V2 + YRS_BTWN_V2INCA, # Calculate Age
         E2_add = sapply(APOE_E2_E3_E4, function(x) sum(as.numeric(unlist(strsplit(as.character(x), ""))) == 2)), # APOE allele count (additive)
         E3_add = sapply(APOE_E2_E3_E4, function(x) sum(as.numeric(unlist(strsplit(as.character(x), ""))) == 3)), 
         E4_add = sapply(APOE_E2_E3_E4, function(x) sum(as.numeric(unlist(strsplit(as.character(x), ""))) == 4)), 
         E2_dom = if_else(is.na(APOE_E2_E3_E4), NA_integer_, if_else(APOE_E2_E3_E4 %in% c(22, 23, 24), 1, 0)),
         E3_dom = if_else(is.na(APOE_E2_E3_E4), NA_integer_, if_else(APOE_E2_E3_E4 %in% c(23, 33, 34), 1, 0)),
         E4_dom = if_else(is.na(APOE_E2_E3_E4), NA_integer_, if_else(APOE_E2_E3_E4 %in% c(24, 34, 44), 1, 0)), 
         GENGROUP = gengrp6, MEAN_AFR = mean_afr, MEAN_EUR = mean_eur, MEAN_AMER = mean_amer) %>%
  select(-c(GENDER, AGE_V2, YRS_BTWN_V2INCA, gengrp6, mean_afr, mean_eur, mean_amer)) %>%
  mutate(across(where(is.character), as.factor)) %>% 
  mutate(APOE_E2_E3_E4 = as.factor(APOE_E2_E3_E4), STRAT = as.factor(STRAT))

# Add labels to factors 
reduced_dat <- reduced_dat %>% 
  mutate(APOE_E2_E3_E4 = factor(APOE_E2_E3_E4, levels = c(22, 23, 24, 33, 34, 44), labels = c("E2_E2", "E2_E3", "E2_E4", "E3_E3", "E3_E4", "E4_E4")), 
         CENTER=factor(CENTER,levels = c("B","C","M","S"),labels=c("Bronx","Chicago", "Miami", "San Diego")), 
         SEX = factor(SEX,levels = c("F","M"),labels=c("Female","Male"))) 

# Creating a missing values plot
vars = c("APOE_E2_E3_E4", 
         "ABETA40", 
         "ABETA42", 
         "PTAU181", 
         "NFLIGHT", 
         "GFAP", 
         "SEX", 
         "AGE",
         "CENTER", 
         "EV1", 
         "MEAN_AFR",
         "GENGROUP")

pdf(file = plot_file, height = 10, width = 10)
gg_miss_upset(reduced_dat[,vars], nsets = length(vars))
dev.off()


# save dataset with missing data 
write.csv(reduced_dat, file = missing_data_file)

# get complete dataset 
# remove censored biomarker data 
complete_dat <- reduced_dat %>% 
  mutate(ABETA40 = if_else(ABETA40_CENSORED == "NC", ABETA40, NA), 
         ABETA42 = if_else(ABETA42_CENSORED == "NC", ABETA42, NA), 
         PTAU181 = if_else(PTAU181_CENSORED == "NC", PTAU181, NA), 
         NFLIGHT = if_else(NFLIGHT_CENSORED == "NC", NFLIGHT, NA), 
         GFAP = if_else(GFAP_CENSORED == "NC", GFAP, NA)) %>% 
  mutate(ABETA40 = if_else(ABETA40 >= 0, ABETA40, NA), 
         ABETA42 = if_else(ABETA42 >= 0, ABETA42, NA), 
         PTAU181 = if_else(PTAU181 >= 0, PTAU181, NA), 
         NFLIGHT = if_else(NFLIGHT >= 0, NFLIGHT, NA), 
         GFAP = if_else(GFAP >= 0, GFAP, NA)) %>% 
  select(-c(ABETA40_CENSORED, ABETA42_CENSORED, GFAP_CENSORED, PTAU181_CENSORED, NFLIGHT_CENSORED))

# want complete APOE genotype and at least one ATN biomarker
complete_dat <- complete_dat %>% 
  filter(!is.na(APOE_E2_E3_E4)) %>% 
  filter(!if_all(all_of(outcome_vars), is.na))


gg_miss_upset(complete_dat[,vars], nsets = length(vars))


# save dataset with "complete" data
write.csv(complete_dat, file = complete_data_file, row.names = FALSE)

