# Creating flow chart for sample selection
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(flowchart)
library(stringr)
library(tidyverse)

# Read in data from csv file 
dat_file = "2024_APOE_ATN_biomarkers/Data/20241001_data.csv"
outliers_removed_dat_file <- "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
outliers_removed <- read.csv(outliers_removed_dat_file) %>% select(SOL_ID, ABETA4240, PTAU181, NFLIGHT, GFAP)

dat <- read.csv(dat_file)
dat <- dat %>% select(-c("ABETA40", "ABETA42", "PTAU181", "NFLIGHT", "GFAP")) %>% 
  mutate(SOL_ID = coalesce(SOL_ID, subject_id, subject_id.x, subject_id.y, SUBJECT_ID.x, SUBJECT_ID.y)) %>%
  left_join(outliers_removed, by = "SOL_ID") %>% 
  mutate(MISSING_APOE_GENOTYPING = if_else(is.na(APOE_E2_E3_E4), TRUE, FALSE), 
         MISSING_CENSORED_BIOMARKER = if_all(all_of(c("ABETA40_CENSORED", "ABETA42_CENSORED", "GFAP_CENSORED", "PTAU181_CENSORED", "NFLIGHT_CENSORED")), ~ is.na(.) | . != "NC"), 
         OUTLIER = if_all(all_of(c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")), ~ is.na(.)) & !MISSING_APOE_GENOTYPING & !MISSING_CENSORED_BIOMARKER, 
         MISSING_GENETIC_BACKGROUND = is.na(gengrp6), 
         MISSING_GENETIC_ANCESTRY = is.na(mean_afr), 
         MISSING_PC = if_any(any_of(c("EV1", "EV2", "EV3", "EV4", "EV5")), ~ is.na(.)), 
         SOL_INCA = if_else(!is.na(WEIGHT_NORM_OVERALL_INCA), TRUE, FALSE), 
         PRIMARY_ANALYSIS = SOL_INCA & !MISSING_APOE_GENOTYPING & !MISSING_CENSORED_BIOMARKER & !OUTLIER, 
         SECONDARY_ANALYSES = PRIMARY_ANALYSIS & !(MISSING_GENETIC_BACKGROUND & MISSING_GENETIC_ANCESTRY) & !MISSING_PC
         ) 


fc <- dat |> 
  as_fc(label = "HCHS/SOL Visit 1", bg_fill = "#F8766D") |> 
  fc_filter(SOL_INCA, label = "SOL/INCA Visit 2", show_exc = TRUE, bg_fill = "#F8766D") |> 
  fc_modify(
    ~ . |> 
      mutate(
        text = ifelse(grepl("Excluded\\n", text), str_glue("{sum(!dat$SOL_INCA)} Excluded from SOL/INCA \n 
                                                           - no follow up at Visit 2 
                                                           - no baseline cognitive data 
                                                           - under age 50"), text)
      )
  ) |> 
  fc_filter(PRIMARY_ANALYSIS, label = "Primary Analysis", show_exc = TRUE, bg_fill = "#619CFF") |> 
  fc_modify(
    ~ . |> 
      mutate(
        text = ifelse(grepl("Excluded\\n", text), str_glue("{sum(dat$SOL_INCA & !dat$PRIMARY_ANALYSIS)} Excluded from Primary Analysis \n
                                                           - {sum(dat$SOL_INCA & dat$MISSING_APOE_GENOTYPING)} missing APOE genotyping 
                                                           - {sum(dat$SOL_INCA & dat$MISSING_CENSORED_BIOMARKER)} missing or censored for all biomarkers
                                                           - {sum(dat$SOL_INCA & dat$OUTLIER)} extreme outliers"), text)
      )
  ) |> 
  fc_filter((SECONDARY_ANALYSES), label = "Secondary Analyses", show_exc = TRUE, bg_fill = "lightblue1") |> 
  fc_modify(
    ~ . |> 
      mutate(
        text = ifelse(grepl("Secondary Analyses", text), str_glue("Secondary Analyses
                                                                  {sum(dat$PRIMARY_ANALYSIS & !dat$MISSING_GENETIC_BACKGROUND & !dat$MISSING_PC)} by Genetic Background Group
                                                                  {sum(dat$PRIMARY_ANALYSIS & !dat$MISSING_GENETIC_ANCESTRY & !dat$MISSING_PC)}  by Genetic Ancestry"), text), 
        text = ifelse(grepl("Excluded\\n", text), str_glue("{sum(dat$PRIMARY_ANALYSIS & !dat$SECONDARY_ANALYSES)} Excluded from Secondary Analyses \n
                                                           - {sum(dat$PRIMARY_ANALYSIS & dat$MISSING_GENETIC_BACKGROUND)} missing genetic background data
                                                           - {sum(dat$PRIMARY_ANALYSIS & dat$MISSING_GENETIC_ANCESTRY)} missing genetic ancestry data
                                                           - {sum(dat$PRIMARY_ANALYSIS & dat$MISSING_PC)} missing genetic PCs"), text)
      )
  ) |> 
  fc_modify(
      ~ . |> 
        mutate(
          bg_fill = ifelse(grepl("Excluded", text), "gray", bg_fill), 
          text_fs = 8
        )
  )


png("2024_APOE_ATN_biomarkers/Results/sample_size_flowchart.png", width=12, height=8, units="in", res=800)
fc |> fc_draw()
dev.off()
