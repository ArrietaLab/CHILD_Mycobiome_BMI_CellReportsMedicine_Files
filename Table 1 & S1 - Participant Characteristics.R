## Participant Characteristics (Table 1 & S1)

# Load packages
library(gtsummary)
library(tidyverse)

# Load metadata for all participants
fungi_data_100 <- read.csv("CHILD 100 Metadata.csv")

# Reorder covariates
fungi_data_100$BF_3m_status <- factor(fungi_data_100$BF_3m_status , levels = c("Zero", "Partial", "Exclusive"))

# Overall participant characteristics table for 100 infants included in this analysis
tbl_summary(fungi_data_100[c("sex", "csec", "Prebirth_abs_oralIV_yn", "Mother_abs_birth_yn", "Child_6moto1Y_abs", "BF_3m_status", "BF_duration_imp", "solids_6m", "site", "hei2010", "AS_bev4", "mom_bmi_best", "Father_BMI")], 
            missing = "no",
            label = c(sex ~ "Infant Sex", csec ~ "Delivery Mode", Prebirth_abs_oralIV_yn ~ "Prenatal Antibiotics", Mother_abs_birth_yn ~ "Intrapartum Antibiotics", BF_3m_status ~ "Breastfeeding Status at 3 Months", solids_6m ~ "Solids Feeding at 6 Months", BF_duration_imp ~ "Breastfeeding Duration (months)", Child_6moto1Y_abs ~ "Infant Antibiotic Exposure (6 months-1 year)", AS_bev4 ~ "Maternal ASB Consumption During Pregnancy", hei2010 ~ "Maternal Healthy Eating Index (2010)", mom_bmi_best ~ "Maternal BMI", Father_BMI ~ "Paternal BMI", site ~ "Study Site"),
            statistic = list(all_continuous() ~ "{mean} ({sd})"),
            digits = all_continuous() ~ 2) %>% modify_header(label ~ "") %>%
  modify_spanning_header(label ~ "**Participant Characteristics**")

# Comparison of BMIz scores by infant sex
tbl_summary(fungi_data_100[c("BMIz_birth", "BMIz_3mo", "BMIz_1yr", "BMIz_3yr", "BMIz_5yr", "sex")], 
            by = "sex",
            missing = "no",
            label = c(BMIz_birth ~ "Birth", BMIz_3mo ~ "3 Months", BMIz_1yr ~ "1 Year", 
                      BMIz_3yr ~ "3 Years", BMIz_5yr ~ "5 Years"),
            statistic = list(all_continuous() ~ "{mean} ({sd})"),
            digits = all_continuous() ~ 2) %>% add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  modify_header(label ~ "**Age**") %>% 
  modify_spanning_header(c("stat_1", "stat_2") ~ "**BMIz Scores in the First 5 Years of Life by Sex**")

# Comparison of covariates evaluated in subsequent analyses by BMI categories (data not shown in manuscript)
tbl_summary(fungi_data_100[c("sex", "csec", "Prebirth_abs_oralIV_yn", "Mother_abs_birth_yn", "Child_6moto1Y_abs", "BF_3m_status", "BF_duration_imp", "solids_6m", "site", "hei2010", "AS_bev4", "mom_bmi_best", "Father_BMI", "BMI_5y_class")], 
            by = "BMI_5y_class",
            missing = "no",
            label = c(sex ~ "Infant Sex", csec ~ "Delivery Mode", Prebirth_abs_oralIV_yn ~ "Prenatal Antibiotics", Mother_abs_birth_yn ~ "Intrapartum Antibiotics", BF_3m_status ~ "Breastfeeding Status at 3 Months", solids_6m ~ "Solids Feeding at 6 Months", BF_duration_imp ~ "Breastfeeding Duration (months)", Child_6moto1Y_abs ~ "Infant Antibiotic Exposure (6 months-1 year)", AS_bev4 ~ "Maternal ASB Consumption During Pregnancy", hei2010 ~ "Maternal Healthy Eating Index (2010)", mom_bmi_best ~ "Maternal BMI", Father_BMI ~ "Paternal BMI", site ~ "Study Site"),
            statistic = list(all_continuous() ~ "{mean} ({sd})"),
            digits = all_continuous() ~ 2) %>% modify_header(label ~ "") %>% add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>% 
  modify_spanning_header(label ~ "**Participant Characteristics by BMI Category at 5 years**")

# Comparison tables repeated for BMI categories at birth, 3 months, 1 year, and 3 years
