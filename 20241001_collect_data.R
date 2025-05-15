require(haven)

# the data that we need:
# APOE genotypes
# ATN biomarkers
# genetic ancestry proportion
# genetic principal components
# Age (at time of ATN biomarkers measures), sex, study center, Hispanic/Latinos background, survey design variables
# (PSU; Strata, SOL-INCA survey weights)

base_folder <- "/Volumes/Sofer Lab/HCHS_SOL/Datasets"
output_folder <- "/Users/tamarsofer/Library/CloudStorage/OneDrive-BethIsraelLaheyHealth/Ongoing_papers/2024_APOE_ATN_biomarkers/Data"

# files with needed data:
apoe_file <- "/Volumes/Sofer Lab/HCHS_SOL/APOE/sol_inca_apoe_inv2.sas7bdat"
ATN_file <- file.path(base_folder, "SOL-INCA/inca2_biomarkers_analysis_2306.sas7bdat")
ancestry_file <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/20240902_global_ancestry_from_2018_local_ancestries/RFMix_global_ancestry_from_LAI.csv"
PC_file <- "/Volumes/Sofer Lab/HCHS_SOL/Ancestry_files/subject_annotation_2017-09-05.csv"

covariates_file <- file.path(base_folder, "sol_cognitive_cov_20210408.csv")
ID_mapping_file <- "/Volumes/Sofer Lab/HCHS_SOL/ID_mapping/hchs_dbgap_id_mapping_200727.sas7bdat"

# read the datasets 
IDs <- read_sas(ID_mapping_file)
IDs$ID <- as.numeric(IDs$ID)

apoe_dat <- read_sas(apoe_file)
dim(apoe_dat) # 12895     6

apoe_dat[1:2,]
# # A tibble: 2 × 6
# SOL_ID    ID       C___3084793_20 C____904973_10 COMBOGT APOE_E2_E3_E4
# <chr>     <chr>    <chr>          <chr>          <chr>   <chr>
#   1 SoL100032 42862192 TT             CC             TT_CC   33
# 2 SoL100100 84506719 TT             CC             TT_CC   33

# turning ID to numeric so we can merge the dataset
# with a dataset that is read from a csv file.
apoe_dat$ID <- as.numeric(apoe_dat$ID)

apoe_dat <- apoe_dat[which(!is.element(apoe_dat$APOE_E2_E3_E4, c(-9, 99))),c("SOL_ID", "ID", "APOE_E2_E3_E4")]


# read ATN data:
ATN <- read_sas(ATN_file)
ATN$ID <- as.numeric(ATN$ID)
nrow(ATN) # 6226
ATN <- merge(ATN, IDs, by = "ID")


dat1 <- merge(apoe_dat, ATN, by = "ID")
nrow(dat1) # 6153
# read ancestry data:

ancestry <- read.csv(ancestry_file)
ancestry$subject_id <- ancestry$SUBJECT_ID
ancestry <- merge(ancestry, IDs, by = "subject_id")
head(ancestry)

dat2 <- merge(dat1, ancestry, by = "ID", all = TRUE)
nrow(dat2)

# read genetic PCs:

pcs <- read.csv(PC_file)
pcs <- pcs[,c("SUBJECT_ID", "EV1", "EV2", "EV3", "EV4", "EV5", "gengrp6")]
pcs$subject_id <- pcs$SUBJECT_ID
pcs <- merge(pcs, IDs, by = "subject_id")


# merge 
dat3 <- merge(dat2, pcs, by ="ID", all = TRUE)
nrow(dat3)


# read covriates and design variables
covars <- read.csv(covariates_file)
covars <- covars[,c("ID", "STRAT", "PSU_ID", "BMI", "GENDER", "AGE_V2", "YRS_BTWN_V2INCA",
                    "WEIGHT_NORM_OVERALL_INCA")]

dat4 <- merge(dat3, covars, by = "ID", all = TRUE)
dim(dat4)

write.csv(dat4, file = file.path(output_folder, "20241001_data.csv"))

