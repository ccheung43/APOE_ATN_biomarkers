# Creating df for Table S2: APOE-ATN associations, stratified by Age or Sex for dominant inheritance mode
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(dplyr)
library(tidyverse)
library(survey)
source("2024_APOE_ATN_biomarkers/Code/format_pvalue.r")
library(openxlsx)
library(readxl)

# Read in data from csv file 
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file) %>%
  mutate(AGE_CAT = case_when(
    AGE < 60 ~ "<60",
    AGE >= 60 & AGE <= 70 ~ "60−70",
    AGE > 70 ~ ">70"
  )) %>% 
  mutate(AGE_CAT = factor(AGE_CAT, levels = c("<60", "60−70", ">70")), 
         SEX = factor(SEX, levels = c("Female", "Male")))


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


APOE_ATN_interaction_by_AGE_SEX <- function(E2_exposure, E4_exposure, modelFormula, num_permutations) { 
  all_tables <- list()
  for(i in 1:nrow(age_sex_groups)) { 
    subset_group <- age_sex_groups[i,]
    print(subset_group)
    svy_obj_subset <- subset(svy_obj, AGE_CAT == subset_group$AGE_CAT & SEX == subset_group$SEX)
    E2 <- vector(mode = "list")
    E4 <- vector(mode = "list")
    for (var in outcome_vars){
      print(var)
      mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = svy_obj_subset, family = gaussian())
      
      E2[[var]] <- extract_one_out(mod, 
                                   exposure = E2_exposure, 
                                   out_name  = var, 
                                   exponentiate = FALSE)
      E4[[var]] <- extract_one_out(mod, 
                                   exposure = E4_exposure, 
                                   out_name  = var, 
                                   exponentiate = FALSE)
      
      # Permutation testing
      e2_permuted_stats <- numeric(num_permutations)
      e4_permuted_stats <- numeric(num_permutations)
      
      for (i in 1:num_permutations) {
        permuted_dat <- dat %>%
          mutate(across(all_of(outcome_vars), ~ ifelse(AGE_CAT == subset_group$AGE_CAT & SEX == subset_group$SEX, sample(.), .))) 
        
        # Fit the model with permuted data
        permuted_svy_obj <- svydesign(id = ~PSU_ID, strata = ~STRAT, weights = ~WEIGHT_NORM_OVERALL_INCA, 
                                      nest = TRUE, data = permuted_dat)
        permuted_svy_obj_subset <- subset(permuted_svy_obj, AGE_CAT == subset_group$AGE_CAT & SEX == subset_group$SEX)
        permuted_mod <- svyglm(formula = as.formula(paste0(var, modelFormula)), design = permuted_svy_obj_subset, family = gaussian())
        
        # Store permuted estimates for E2 and E4 exposures
        e2_permuted_stats[i] <- summary(permuted_mod)$coef[E2_exposure, "Estimate"]
        e4_permuted_stats[i] <- summary(permuted_mod)$coef[E4_exposure, "Estimate"]
      }
      
      # Calculate permutation-based p-values
      if (E2[[var]]$est < 0) {
        E2[[var]]$perm_pval <- mean(e2_permuted_stats <= E2[[var]]$est)
      } else {
        E2[[var]]$perm_pval <- mean(e2_permuted_stats >= E2[[var]]$est)
      }
      
      if (E4[[var]]$est < 0) {
        E4[[var]]$perm_pval <- mean(e4_permuted_stats <= E4[[var]]$est)
      } else {
        E4[[var]]$perm_pval <- mean(e4_permuted_stats >= E4[[var]]$est)
      }
    } 
    
    E2_tab <- do.call(rbind, E2) %>% mutate(allele = "E2")
    E4_tab <- do.call(rbind, E4) %>% mutate(allele = "E4")
    
    subset_table <- rbind(E2_tab, E4_tab) %>% 
      mutate(biomarker = factor(outcome, levels = c("ABETA4240", "PTAU181","NFLIGHT", "GFAP"))) %>% 
      arrange(biomarker) %>% 
      select(biomarker, allele, est, CI, perm_pval) %>% 
      mutate(paper_pval = format_pvalue(perm_pval))
    
    group_name = paste0(subset_group$SEX, "_", subset_group$AGE_CAT)
    all_tables[[group_name]] <- subset_table
  } 
  return(all_tables)
} 

outcome_vars <- c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")
age_sex_groups <- expand.grid(AGE_CAT = levels(dat$AGE_CAT), SEX = levels(dat$SEX))  

# Additive inheritance mode, stratified by AGE and SEX
e2_exposure = "E2_dom"
e4_exposure = "E4_dom"
model_formula = "~E2_dom+E4_dom+CENTER"
dom_tbl_by_age_sex <- APOE_ATN_interaction_by_AGE_SEX(e2_exposure, e4_exposure, model_formula, 10000)
association_file = "2024_APOE_ATN_biomarkers/Data/20250331_age_sex_stratified_dom_table.xlsx" 
# Write to xlsx sheet with a page per age group
write.xlsx(dom_tbl_by_age_sex, file = association_file, overwrite = TRUE) 



# write function to format estimate as est[CI]
format_est_ci <- function(est, CI) {
  # Extract lower and upper bounds from CI
  bounds <- as.numeric(unlist(strsplit(gsub("[()]", "", CI), ",")))
  # Combine into single string
  return(sprintf("%s[%s,%s]", est, bounds[1], bounds[2]))
}

dat_file <- "2024_APOE_ATN_biomarkers/Data/20250331_age_sex_stratified_dom_table.xlsx"
sheet_names <- excel_sheets(dat_file)
associations_tables <- lapply(sheet_names, function(sheet) read_excel(dat_file, sheet = sheet))
names(associations_tables) <- sheet_names 


# stratify by age first, then sex
tbls2 <- data.frame("ATN biomarker" = c("Aꞵ42/40", "", "pTau-181", "", "NfL", "", "GFAP", ""), 
                    "APOE allele" = c("ε2", "ε4", "ε2", "ε4", "ε2", "ε4", "ε2", "ε4"))

sheet_names <- c("Female_<60", "Male_<60", "Female_60−70", "Male_60−70", "Female_>70", "Male_>70")
for(group in sheet_names) { 
  group_df <- associations_tables[[group]]
  group_df <- group_df  %>% 
    mutate("est[CI]" = mapply(format_est_ci, est, CI)) %>% 
    mutate("p-value" = paper_pval) %>% 
    select(c("est[CI]", "p-value"))
  tbls2 <- cbind(tbls2, group_df)
} 
row1 <- c("", "", "<60", "", "", "","60-70", "", "", "", ">70", "", "", "")
row2 <- c("", "", "Female", "", "Male", "", "Female", "", "Male", "", "Female", "", "Male", "")
row3 <- c("ATN biomarker", "APOE allele", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value")
tbls2 <- rbind(row1, row2, row3, tbls2)


tbl_file = "2024_APOE_ATN_biomarkers/Results/20250331_table_s2.xlsx"
write.xlsx(tbls2, file = tbl_file, overwrite = TRUE, colNames=FALSE) 
