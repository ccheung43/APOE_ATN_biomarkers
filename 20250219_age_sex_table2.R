# Create Table 3: Associations between ATN biomarkers and APOE allele by age and sex
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(openxlsx)
library(readxl)

# write function to format estimate as est[CI]
format_est_ci <- function(est, CI) {
  # Extract lower and upper bounds from CI
  bounds <- as.numeric(unlist(strsplit(gsub("[()]", "", CI), ",")))
  # Combine into single string
  return(sprintf("%s[%s,%s]", est, bounds[1], bounds[2]))
}

dat_file <- "2024_APOE_ATN_biomarkers/Data/20250219_age_sex_stratified_add_table.xlsx" 
sheet_names <- excel_sheets(dat_file)
associations_tables <- lapply(sheet_names, function(sheet) read_excel(dat_file, sheet = sheet))
names(associations_tables) <- sheet_names 


# stratify by sex first, then age
tbl3a <- data.frame("ATN biomarker" = c("Aꞵ42/40", "", "pTau-181", "", "NfL", "", "GFAP", ""), 
                   "APOE allele" = c("ε2", "ε4", "ε2", "ε4", "ε2", "ε4", "ε2", "ε4"))
for(group in sheet_names) { 
  group_df <- associations_tables[[group]]
  group_df <- group_df  %>% 
    mutate("est[CI]" = mapply(format_est_ci, est, CI)) %>% 
    mutate("p-value" = paper_pval) %>% 
    select(c("est[CI]", "p-value"))
  tbl3a <- cbind(tbl3a, group_df)
} 
row1 <- c("", "", "Female", "", "", "","", "", "Male", "", "", "", "", "")
row2 <- c("", "", "<60", "", "60-70", "", ">70", "", "<60", "", "60-70", "", ">70", "")
row3 <- c("ATN biomarker", "APOE allele", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value")
tbl3a <- rbind(row1, row2, row3, tbl3a)


tbl3_file = "2024_APOE_ATN_biomarkers/Results/20250219_table3_age_sex.xlsx"
write.xlsx(tbl3a, file = tbl3_file, overwrite = TRUE, colNames=FALSE)  

# stratify by age first, then sex
tbl3b <- data.frame("ATN biomarker" = c("Aꞵ42/40", "", "pTau-181", "", "NfL", "", "GFAP", ""), 
                    "APOE allele" = c("ε2", "ε4", "ε2", "ε4", "ε2", "ε4", "ε2", "ε4"))

sheet_names <- c("Female_<60", "Male_<60", "Female_60−70", "Male_60−70", "Female_>70", "Male_>70")
for(group in sheet_names) { 
  group_df <- associations_tables[[group]]
  group_df <- group_df  %>% 
    mutate("est[CI]" = mapply(format_est_ci, est, CI)) %>% 
    mutate("p-value" = paper_pval) %>% 
    select(c("est[CI]", "p-value"))
  tbl3b <- cbind(tbl3b, group_df)
} 
row1 <- c("", "", "<60", "", "", "","60-70", "", "", "", ">70", "", "", "")
row2 <- c("", "", "Female", "", "Male", "", "Female", "", "Male", "", "Female", "", "Male", "")
row3 <- c("ATN biomarker", "APOE allele", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value", "est[CI]", "p-value")
tbl3b <- rbind(row1, row2, row3, tbl3b)


tbl3_file = "2024_APOE_ATN_biomarkers/Results/20250311_table3_age_sex.xlsx"
write.xlsx(tbl3b, file = tbl3_file, overwrite = TRUE, colNames=FALSE) 
