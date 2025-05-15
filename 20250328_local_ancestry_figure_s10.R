# Creating a figure for Figure S10: Interaction of local genetic ancestry counts and APOE ε2 allele on ATN biomarkers. 
# by testing linear combination of regression coefficients for proportion genetic ancestry
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(tidyverse)
library(readxl)
library(stringr) 
library(Cairo)

# file names
dat_file = "2024_APOE_ATN_biomarkers/Data/20250124_local_genetic_ancestry_linear_combinations.xlsx"
fig_file = "2024_APOE_ATN_biomarkers/Results/20250328_figure_s10.png"
# Read in data from xlsx file 
df <- read_excel(dat_file)

df <- df %>% 
  mutate(est = signif(as.numeric(combined_est), 3)) %>% 
  mutate(CI_Lower = signif(est - 2*as.numeric(combined_est_se)), 
         CI_Upper = signif(est + 2*as.numeric(combined_est_se)), 
         significance = case_when(
           pval < 0.001 ~ "***",
           pval < 0.01 ~ "**",
           pval < 0.05 ~ "*",
           TRUE ~ ""), 
         ancestry = str_to_title(ancestry)) %>% 
  mutate(biomarker = factor(biomarker, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")))

fig <- df %>% 
  filter(allele == "E2") %>% 
  mutate(counts = as.character(as.numeric(proportion) * 2)) %>% 
  ggplot(aes(x = counts, y = est, group = ancestry, color = ancestry)) +
  geom_line(linewidth = 1.2, alpha = 0.5) + 
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = ancestry), alpha = 0.2, color = NA, show.legend = FALSE) + 
  geom_point(size = 3) + 
  facet_wrap(~biomarker, nrow = 2, scale = "free_y",
             labeller = labeller(biomarker = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP"), 
             )) +
  labs(x = "Local Genetic Ancestry Counts", y = "Estimate and 95% CI", 
       color = "Genetic Ancestry") + 
  theme_bw(base_size = 16)

#CairoPDF(fig_file, width = 16, height = 6)
CairoPNG(fig_file, width = 10, height = 6, units = "in", res = 300)
print(fig)
dev.off()
