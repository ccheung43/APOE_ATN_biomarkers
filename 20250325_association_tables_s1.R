# Creating Tables S1: Associations between ATN biomarkers and APOE allele
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(dplyr)
library(tidyverse)
library(survey)
source("2024_APOE_ATN_biomarkers/Code/format_pvalue.r")
library(openxlsx)

# Read in data from csv file 
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file)
association_file = "2024_APOE_ATN_biomarkers/Data/20250325_associations_tables.xlsx" 
s1_file = "2024_APOE_ATN_biomarkers/Results/20250325_s1.xlsx"

# survey object
svy_obj <- svydesign(id = ~PSU_ID,
                     strata = ~STRAT,
                     weights = ~WEIGHT_NORM_OVERALL_INCA,
                     nest = TRUE,
                     data = dat)

# from Gendered_indices_and_insomnia 
extract_one_out <- function(mod, exposure, out_name, exponentiate = FALSE, sigfigs = 2){
  est <- summary(mod)$coef
  confint <- confint(mod)
  
  ind <- grep(exposure, rownames(est))
  
  if (exponentiate){
    df_out <- data.frame(outcome = out_name, 
                         est = signif(exp(est[ind, "Estimate"]),sigfigs), 
                         CI = paste0("(", signif(exp(confint[ind, 1]), sigfigs), ",",
                                     signif(exp(confint[ind, 2]), sigfigs), ")"),
                         pval = formatC(est[ind, "Pr(>|t|)"], digits = sigfigs, format = "E"),
                         paper_pval= format_pvalue(est[ind, "Pr(>|t|)"]))
  } else{
    df_out <- data.frame(outcome = out_name, 
                         est = signif(est[ind, "Estimate"],sigfigs), 
                         CI = paste0("(", signif(confint[ind, 1], sigfigs), ",",
                                     signif(confint[ind, 2], sigfigs), ")"),
                         pval = formatC(est[ind, "Pr(>|t|)"], digits = sigfigs, format = "E"),
                         paper_pval = format_pvalue(est[ind, "Pr(>|t|)"]))
  }
  
  return(df_out)
}
# Write a function to extract the svyglm estimates, CIs, and pvalues for the 2 APOE allele exposures and 4 ATN biomarker outcomes
extract_associations_table <- function(E2_exposure, E4_exposure, modelFormula) {
  outcome_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
  E2 <- vector(mode = "list")
  E4 <- vector(mode = "list")
  for (var in outcome_vars){
    mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = svy_obj, family = gaussian())
    E2[[var]] <- extract_one_out(mod, 
                                 exposure = E2_exposure, 
                                 out_name  = var, 
                                 exponentiate = FALSE)
    E4[[var]] <- extract_one_out(mod, 
                                 exposure = E4_exposure, 
                                 out_name  = var, 
                                 exponentiate = FALSE) 
  }
  
  E2_tab <- do.call(rbind, E2) %>% mutate(allele = "E2")
  E4_tab <- do.call(rbind, E4) %>% mutate(allele = "E4")
  associations_table <- rbind(E2_tab, E4_tab) %>% 
    mutate(biomarker = factor(outcome, levels = c("ABETA4240", "PTAU181","NFLIGHT", "GFAP"))) %>% 
    arrange(biomarker) %>% 
    select(biomarker, allele, est, CI, pval, paper_pval)
}

# Create a list to store each table
associations_tables <- list()

# Dominant inheritance mode (unadjusted table)
E2_exposure = "E2_dom"
E4_exposure = "E4_dom"
modelFormula = "~E2_dom+E4_dom"
associations_tables[["DominantModel0"]] <- extract_associations_table(E2_exposure, E4_exposure, modelFormula)

# Dominant inheritance mode (patrially adjusted table)
E2_exposure = "E2_dom"
E4_exposure = "E4_dom"
modelFormula =  "~E2_dom+E4_dom+SEX+AGE+CENTER"
associations_tables[["DominantModel1"]] <- extract_associations_table(E2_exposure, E4_exposure, modelFormula)

# Dominant inheritance mode (adjusted table)
E2_exposure = "E2_dom"
E4_exposure = "E4_dom"
modelFormula =  "~E2_dom+E4_dom+SEX+AGE+CENTER+EV1+EV2+EV3+EV4+EV5"
associations_tables[["DominantModel2"]] <- extract_associations_table(E2_exposure, E4_exposure, modelFormula)


# Additive inheritance mode (unadjusted table)
E2_exposure = "E2_add"
E4_exposure = "E4_add"
modelFormula = "~E2_add+E4_add"
associations_tables[["AdditiveModel0"]] <- extract_associations_table(E2_exposure, E4_exposure, modelFormula)

# Additive inheritance mode (partially adjusted table)
E2_exposure = "E2_add"
E4_exposure = "E4_add"
modelFormula =  "~E2_add+E4_add+SEX+AGE+CENTER"
associations_tables[["AdditiveModel1"]] <- extract_associations_table(E2_exposure, E4_exposure, modelFormula)

# Additive inheritance mode (adjusted table)
E2_exposure = "E2_add"
E4_exposure = "E4_add"
modelFormula =  "~E2_add+E4_add+SEX+AGE+CENTER+EV1+EV2+EV3+EV4+EV5"
associations_tables[["AdditiveModel2"]] <- extract_associations_table(E2_exposure, E4_exposure, modelFormula)


# Write each table to a separate page of a excel sheet 
wb <- createWorkbook()
for (i in seq_along(associations_tables)) {
  addWorksheet(wb, sheetName = names(associations_tables)[i])  # Use list names as sheet names
  writeData(wb, sheet = i, associations_tables[[i]])
}
saveWorkbook(wb, association_file, overwrite = TRUE)


# Format Table S! and write to xlsx 

# write function to format table2 estimate as est[CI]
format_est_ci <- function(est, CI) {
  # Extract lower and upper bounds from CI
  bounds <- as.numeric(unlist(strsplit(gsub("[()]", "", CI), ",")))
  # Combine into single string
  return(sprintf("%s[%s,%s]", est, bounds[1], bounds[2]))
}

s1 <- data.frame("ATN biomarker" = c("Aꞵ42/40", "", "pTau-181", "", "NfL", "", "GFAP", ""), 
                    "APOE allele" = c("ε2", "ε4", "ε2", "ε4", "ε2", "ε4", "ε2", "ε4"))

sheet_names <- c("AdditiveModel0", "DominantModel0", "AdditiveModel1","DominantModel1", "AdditiveModel2", "DominantModel2")
for(mod in sheet_names) { 
  mod_df <- associations_tables[[mod]]
  mod_df <- mod_df  %>% 
    mutate("est[CI]" = mapply(format_est_ci, est, CI)) %>% 
    mutate("p-value" = paper_pval) %>% 
    select(c("est[CI]", "p-value"))
  s1 <- cbind(s1, mod_df)
} 


row1 <- c("", "", "Model 0", "", "", "", "Model 1", "", "", "", "Model2", "", "", "")
row2 <- c("", "", "Additive", "", "Dominant", "", "Additive", "", "Dominant", "", "Additive", "", "Dominant", "")
row3 <- c("ATN biomarker", "APOE allele", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value") 
s1 <- rbind(row1, row2, row3, s1)


write.xlsx(s1, file = s1_file, overwrite = TRUE, colNames=FALSE) 


