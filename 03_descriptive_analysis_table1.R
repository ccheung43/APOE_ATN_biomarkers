# Creating df for Table 1 (with outliers removed): Demographics, genetics, and biomarkers of HCHS/SOL by genetic background groups
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(dplyr)
library(tidyverse)
library(tableone)
library(survey)
library(forcats)
library(openxlsx)

# Read in files
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file)
table_file = "2024_APOE_ATN_biomarkers/Results/20250121_table1_outliers_removed.xlsx"


dat <- dat %>%
  mutate(AGE_CAT = case_when(
    AGE < 60 ~ "<60",
    AGE >= 60 & AGE <= 70 ~ "60−70",
    AGE > 70 ~ ">70"
  )) %>% 
  mutate(AGE_CAT = factor(AGE_CAT, levels = c("<60", "60−70", ">70"))) %>% 
  mutate(APOE_E2_E3_E4 = fct_recode(APOE_E2_E3_E4,
                                    "ε2/ε2" = "E2_E2",
                                    "ε2/ε3" = "E2_E3",
                                    "ε2/ε4" = "E2_E4", 
                                    "ε3/ε3" = "E3_E3",
                                    "ε3/ε4" = "E3_E4",
                                    "ε4/ε4" = "E4_E4"))
vars <-  c("SEX", 
           "AGE", 
           "AGE_CAT", 
           "APOE_E2_E3_E4", 
           "ABETA40", 
           "ABETA42", 
           "PTAU181", 
           "NFLIGHT", 
           "GFAP", 
           "MEAN_AFR", 
           "MEAN_EUR", 
           "MEAN_AMER")

survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_INCA , data=dat)

tbl1_unweighted <- CreateTableOne(vars = vars, data = dat, strata = "GENGROUP", addOverall = TRUE, test = FALSE) %>%
  print(showAllLevels = TRUE, varLabels = TRUE, digits = 3)

tbl1_weighted <- svyCreateTableOne(vars = vars, data = survey_obj, strata = "GENGROUP", addOverall = TRUE, test = FALSE) %>%
  print(showAllLevels = TRUE, varLabels = TRUE, digits = 3)

# table for alleles E2, E3, and E4
dat_alleles <- dat %>% 
  pivot_longer(c(E2_add, E3_add, E4_add), names_to = "APOE_allele", values_to = "APOE_allele_count") %>% 
  uncount(APOE_allele_count) %>% 
  mutate(APOE_allele = as.factor(APOE_allele)) %>% 
  mutate(APOE_allele = fct_recode(APOE_allele,"ε2"="E2_add", "ε3"="E3_add", "ε4"="E4_add"))

survey_obj <- svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_INCA , data=dat_alleles)

alleles_unweighted <- CreateTableOne(vars = "APOE_allele", data = dat_alleles, strata = "GENGROUP", addOverall = TRUE, test = FALSE) %>%
  print(showAllLevels = TRUE, varLabels = TRUE, digits = 3)

alleles_weighted <- svyCreateTableOne(vars = "APOE_allele", data = survey_obj, strata = "GENGROUP", addOverall = TRUE, test = FALSE) %>%
  print(showAllLevels = TRUE, varLabels = TRUE, digits = 3)

tbl1_unweighted <- rbind(tbl1_unweighted[1:13,], alleles_unweighted[2:4,], tbl1_unweighted[14:21,])
tbl1_weighted <- rbind(tbl1_weighted[1:13,], alleles_weighted[2:4,], tbl1_weighted[14:21,])

tbl1_comb <- tbl1_unweighted
col_inds_to_update <- 2:8 # all columns except "level"
row_inds_to_update <- c(1:3, 5:13) # rows corresponding to N, Sex, Age (categorical), APOE genotype

# create function to remove any leading/trailing spaces around the parentheses
fix_spacing <- function(x) {
  gsub("\\s*\\(", " (", gsub("\\s*\\)", ") ", gsub("\\s*\\(\\s*", " (", trimws(x)))) 
}

# update tbl1_comb with the percentages from the weighted table
for (i in col_inds_to_update){
  counts <- sapply(tbl1_unweighted[,i], function(x){
    strsplit(x, split = "(", fixed = TRUE)[[1]][1]
  })
  percent <- sapply(tbl1_weighted[,i], function(x){
    paste0("(",  strsplit(x, split = "(", fixed = TRUE)[[1]][2])
  })
  percent[which(names(percent) == "n")] = ""
  tbl1_comb[,i] <- paste0(counts, percent)
  tbl1_comb[,i] <- fix_spacing(tbl1_comb[, i])
}
rownames(tbl1_comb) = c("N", 
                        "Sex (%)", "", 
                        "Age (mean (SD))", 
                        "Age (%)", "", "", 
                        "APOE genotype (%)",  "", "", "", "", "", 
                        "APOE allele (%)", "", "",
                        "Aꞵ40 (mean (SD))", "Aꞵ42 (mean (SD))", "pTau-181 (mean (SD))", "NfL (mean (SD))", "GFAP (mean (SD))", 
                        "African (mean (SD))", "European (mean(SD))", "Amerindian (mean (SD))")

tbl1 = cbind(rownames(tbl1_comb), tbl1_comb)

write.xlsx(tbl1, file = table_file, overwrite = TRUE)


t1flex(tbl1)

