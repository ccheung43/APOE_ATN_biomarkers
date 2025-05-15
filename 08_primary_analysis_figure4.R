# Creating figure for Table 2 to show a forest plot of each model for each biomarker and model
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(readxl)
library(stringr) 
library(Cairo)

# file names
dat_file = "2024_APOE_ATN_biomarkers/Data/20250325_associations_tables.xlsx" 
fig_file = "2024_APOE_ATN_biomarkers/Results/20250313_figure4.png"

# Read in data from xlsx file 
sheet_names <- excel_sheets(dat_file)
all_sheets <- lapply(sheet_names, function(sheet) read_excel(dat_file, sheet = sheet))
names(all_sheets) <- sheet_names 


df <- data.frame()
for(i in 1:length(all_sheets)) { 
  sheet <- as.data.frame(all_sheets[[i]]) %>% mutate(mod = names(all_sheets)[i])
  df <- rbind(df, sheet)
}

fig <- df %>% 
  filter(str_detect(mod, "Additive")) %>% 
  mutate(mod = case_when(
    mod == "AdditiveModel0" ~ "Model 0", 
    mod == "AdditiveModel1" ~ "Model 1", 
    mod == "AdditiveModel2" ~ "Model 2")) %>% 
  mutate(allele = if_else(allele == "E2", paste0("\u03B5", "2"), paste0("\u03B5", "4"))) %>% 
  mutate(pval = as.numeric(pval)) %>% 
  mutate(CI = str_remove_all(CI, "[()]"),              
         CI_Lower = as.numeric(str_extract(CI, "^[^,]+")), 
         CI_Upper = as.numeric(str_extract(CI, "[^,]+$")), 
         significance = case_when(
           pval < 0.001 ~ "***",
           pval < 0.01 ~ "**",
           pval < 0.05 ~ "*",
           TRUE ~ "")) %>% 
  mutate(biomarker = factor(biomarker, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) %>% 
  ggplot(aes(x = allele, y = est, color = mod)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), linewidth = 1, width = 0.2, alpha = 0.5, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = significance, group = mod), size = 4, color = "black", position = position_dodge(width = 0.5)) + 
  labs(y = "Estimate and 95% CI", x = "APOE allele", color = "Model", shape = "APOE allele") +
  facet_wrap(~biomarker, nrow = 2, scales = "free_y", 
             labeller = labeller(biomarker = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP")))+ 
  #scale_shape_discrete(labels = c(paste0("\u03B5", "2"), paste0("\u03B5", "4"))) + 
  theme_bw(base_size = 14)

#CairoPDF(fig_file, width = 16, height = 4)  # Specify size in inches for PDF
CairoPNG(fig_file, width = 8, height = 4, units = "in", res = 300)
print(fig)
dev.off()
