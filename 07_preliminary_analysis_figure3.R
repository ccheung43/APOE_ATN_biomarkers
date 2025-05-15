# Creating Figure 3 to show ATN biomarker levels by APOE genotype
# Assuming working directory is OneDrive - Beth Israel Lahey Health

# Load required libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Cairo)
library(patchwork)

# Read in data from csv file 
dat_file = "2024_APOE_ATN_biomarkers/Data/20241031_data_outliers_removed.csv"
dat <- read.csv(dat_file)
fig_file = "2024_APOE_ATN_biomarkers/Results/20250313_figure3.png" #change to .pdf if you want a pdf file

fig_top <- dat %>%
  mutate(APOE_E2_E3_E4 = fct_recode(APOE_E2_E3_E4,
                                    "ε2/ε2" = "E2_E2",
                                    "ε2/ε3" = "E2_E3",
                                    "ε2/ε4" = "E2_E4", 
                                    "ε3/ε3" = "E3_E3",
                                    "ε3/ε4" = "E3_E4",
                                    "ε4/ε4" = "E4_E4")) %>% 
  pivot_longer(c(ABETA4240, PTAU181, NFLIGHT, GFAP), 
               names_to = "BIOMARKER", values_to = "CONCENTRATION") %>%
  mutate(BIOMARKER = factor(BIOMARKER, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP"))) %>%
  group_by(APOE_E2_E3_E4, BIOMARKER) %>%
  summarise(mean_concentration = mean(CONCENTRATION, na.rm = TRUE),
            sd_concentration = sd(CONCENTRATION, na.rm = TRUE),
            .groups = 'drop') %>% 
  ggplot(aes(x = as.factor(APOE_E2_E3_E4), y = mean_concentration, color = APOE_E2_E3_E4)) +
  geom_point(stat = "identity", position = position_dodge(0.5), size = 3) + 
  geom_errorbar(aes(ymin = mean_concentration - sd_concentration, 
                    ymax = mean_concentration + sd_concentration),
                linewidth = 1, width = 0.2, position = position_dodge(0.5), alpha = 0.5) +
  labs(y = "", x = "APOE Genotype", color = "Genotype") +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_blank())  + 
  facet_wrap(~BIOMARKER, scales = "free_y", nrow = 1, 
             labeller = labeller(BIOMARKER = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40*"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP")))

fig_bottom <- dat %>% 
  pivot_longer(c(ABETA4240, PTAU181, NFLIGHT, GFAP), 
               names_to = "BIOMARKER", values_to = "CONCENTRATION") %>%
  mutate(BIOMARKER = factor(BIOMARKER, levels = c("ABETA4240", "PTAU181", "NFLIGHT", "GFAP")), 
         E4_add = as.factor(E4_add)) %>%
  group_by(E4_add, BIOMARKER) %>%
  summarise(mean_concentration = mean(CONCENTRATION, na.rm = TRUE),
            sd_concentration = sd(CONCENTRATION, na.rm = TRUE),
            .groups = 'drop') %>% 
  ggplot(aes(x = E4_add, y = mean_concentration, color = E4_add)) +
  geom_point(stat = "identity", position = position_dodge(0.5), size = 3) + 
  geom_errorbar(aes(ymin = mean_concentration - sd_concentration, 
                    ymax = mean_concentration + sd_concentration),
                linewidth = 1, width = 0.2, position = position_dodge(0.5), alpha = 0.5) +
  labs(y = "                                  Mean Concentration (pg/ml)", 
       x = "APOE ε4 Count", color = "Count") +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_blank())  +
  facet_wrap(~BIOMARKER, scales = "free_y", nrow = 1, 
             labeller = labeller(BIOMARKER = c(
               "ABETA4240" = paste0("A", "\u03B2","42/40*"), 
               "PTAU181" = "pTau-181", 
               "NFLIGHT" = "NfL", 
               "GFAP" = "GFAP")))
fig <- fig_top / fig_bottom

#CairoPDF(fig_file, width = 10, height = 5)  # Specify size in inches for PDF
CairoPNG(fig_file, width = 8, height = 4, units = "in", res = 300)
print(fig)
dev.off()



