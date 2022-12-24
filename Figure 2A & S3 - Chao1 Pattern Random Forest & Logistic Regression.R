## Chao1 Pattern Random Forest & Logistic Regression (Figures 2A & S3)

# Load packages
library(tidyverse)
library(phyloseq)

# Load phyloseq object
load("ps_fungi_clean")

# Extract sample data from phyloseq object to use for analysis 
fungi_data <- as(sample_data(ps_fungi_clean), "data.frame")

# Create filtered dataframe by one time point only for analysis
fungi_data_3 <- fungi_data %>%
  filter(Sample_time == "3")

# Remove NA and unchanged values for Chao1 pattern from dataframe
sdf2 <- fungi_data_3 %>% filter(Chao1_dif_cat == c("Increase","Decrease"))

# Prepare data for random forest analysis. Before running model, we must choose the level of our covariates that we want to use as our baseline.
sdf2$Chao1_dif_cat <- as.factor(sdf2$Chao1_dif_cat)
sdf2$Chao1_dif_cat <- relevel(sdf2$Chao1_dif_cat, ref = "Increase")

sdf2$sex <- as.factor(sdf2$sex)
sdf2$sex <- relevel(sdf2$sex, ref = "Male")

sdf2$csec <- as.factor(sdf2$csec)
sdf2$csec <- relevel(sdf2$csec, ref = "Vaginal")

sdf2$BF_3m_status_new <- as.factor(sdf2$BF_3m_status_new)
sdf2$BF_3m_status_new <- relevel(sdf2$BF_3m_status_new, ref = "Zero/Partial")

sdf2$solids_6m <- as.factor(sdf2$solids_6m)
sdf2$solids_6m <- relevel(sdf2$solids_6m, ref = "No")

sdf2$Prebirth_abs_oralIV_yn <- as.factor(sdf2$Prebirth_abs_oralIV_yn)
sdf2$Prebirth_abs_oralIV_yn <- relevel(sdf2$Prebirth_abs_oralIV_yn, ref = "No")

sdf2$Mother_abs_birth_yn <- as.factor(sdf2$Mother_abs_birth_yn)
sdf2$Mother_abs_birth_yn <- relevel(sdf2$Mother_abs_birth_yn, ref = "No")

sdf2$Child_6moto1Y_abs <- as.factor(sdf2$Child_6moto1Y_abs)
sdf2$AS_bev4 <- as.factor(sdf2$AS_bev4)

sdf2$Site_Vancouver <- as.factor(sdf2$Site_Vancouver)
sdf2$Site_Vancouver <- relevel(sdf2$Site_Vancouver, ref = "Vancouver")

# Create new dataframe only containing covariates to be used in subsequent analyses
x <- sdf2[, c("Chao1_dif_cat", "sex", "csec", "BF_3m_status_new", "Site_Vancouver", "solids_6m", "BF_duration_imp", "mom_bmi_best", "Father_BMI", "Mother_abs_birth_yn", "Prebirth_abs_oralIV_yn", "Child_6moto1Y_abs", "hei2010", "AS_bev4","BMIz_3mo", "BMIz_1yr", "Shannon_bacteria_3m", "PCoA1_bacteria_3m", "Shannon_bacteria_12m", "PCoA1_bacteria_12m")]

# Determine number of NA's per covariate
summary(x)

# Load packages required for imputation
library(mice)
library(lattice)

# Set seed for reproducible results
set.seed(123)

# Check the missingness pattern for the dataset
m <- md.pattern(x)
m

# Impute data for paternal BMI and BMIz at 3 months and 1 year to prevent decreased sample size for random forest
imp <- mice(x, seed = 123)
imp
imp$imp

c.broad <- complete(imp, "broad")
c.broad <- c.broad %>% 
  mutate(Father_BMI_imp = (Father_BMI.1 + Father_BMI.2 + Father_BMI.3 + Father_BMI.4 + Father_BMI.5)/5) %>% 
  mutate(BMIz_3mo_imp = (BMIz_3mo.1 + BMIz_3mo.2 + BMIz_3mo.3 + BMIz_3mo.4 + BMIz_3mo.5)/5) %>% 
  mutate(BMIz_1yr_imp = (BMIz_1yr.1 + BMIz_1yr.2 + BMIz_1yr.3 + BMIz_1yr.4 + BMIz_1yr.5)/5)

# Check imputation
print(c.broad$Father_BMI_imp)
print(c.broad$BMIz_3yr_imp)
print(c.broad$BMIz_5yr_imp)

# Replace values with imputed data
x$Father_BMI <- c.broad$Father_BMI_imp
x$BMIz_3mo <- c.broad$BMIz_3mo_imp
x$BMIz_1yr <- c.broad$BMIz_1yr_imp

# Remove NAs from newly imputed dataframe
x <- x %>% na.omit

# Check dataframe for sample size and value distribution
summary(x)

# Load packages required for random forest analysis
library(randomForest)
library(caret)

# Set seed for reproducible results
set.seed(151) 

# Train dataset for random forest
fit_control <- trainControl(method = "cv", number = 10 )

# Generate random forest to determine factors that are most predictive/strongly associated with Chao1 pattern
RF_state_classify <- randomForest(Chao1_dif_cat ~ . , data=x, importance=TRUE, proximity=TRUE, nperm=1000, ntree=500, trControl=fit_control)
RF_state_classify

# Determine feature importance and sort by mean decreasing Gini index
RF_state_classify_imp <- as.data.frame(RF_state_classify$importance)
RF_state_classify_imp$features <- rownames(RF_state_classify_imp)
RF_state_classify_imp_sorted <- arrange(RF_state_classify_imp  , desc(MeanDecreaseGini))
m <- as.data.frame(as.matrix(RF_state_classify_imp_sorted))

# Create rounding function for rounding numeric values in dataframe
round_df <- function(x, digits) {
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

# Round mean decreasing Gini index
m$MeanDecreaseGini <- as.numeric(m$MeanDecreaseGini)
m <- round_df(m, 2)

# Assign features to categories for colour coding figure
m1 <- m %>% mutate(Factor = case_when(features %in% c("sex", "csec", "Site_Vancouver", "Child_6moto1Y_abs", "Mother_abs_birth_yn", "age_stool") ~ "Early Life Factors", features %in% c("BF_3m_status_new", "BF_duration_imp", "solids_6m") ~ "Infant Nutrition", features %in% c("Shannon_bacteria_3m","Shannon_bacteria_12m","PCoA1_bacteria_3m","PCoA1_bacteria_12m") ~ "Infant Microbiome", features %in% c("mom_bmi_best", "hei2010", "Prebirth_abs_oralIV_yn", "AS_bev4") ~ "Maternal Factors", features %in% c("BMIz_3mo", "BMIz_1yr") ~ "Infant BMIz", features %in% c("Father_BMI") ~ "Paternal Factors"))

# Plot random forest results by mean decreasing Gini index (Figure 2A)
random_forest_fig <- ggplot(m1, aes(reorder(features, MeanDecreaseGini), MeanDecreaseGini, fill=Factor))+ 
  geom_bar(stat="identity",  color = "black", width=0.6)+
  labs(colour = NULL)+ labs(x="", y ="Mean Decreasing Gini Index")+
  coord_flip()+
  scale_fill_manual(values=c("purple4","mediumorchid4","deepskyblue3","deepskyblue4","#35b779","olivedrab3"), name = "Factor")+
  theme_bw()+
  ggtitle("Predictors of Fungal Richness Pattern")+
  theme(plot.title = element_text(size=10, hjust = 0.5, face = "bold"),
        panel.background = element_blank(),
        axis.title=element_text(face = "bold", size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9),
        legend.text.align = 0)+
  scale_x_discrete(labels=c("hei2010" = "Maternal Healthy Eating Index", "BMIz_3mo" = "Infant BMIz (3 months)", "mom_bmi_best" = "Maternal BMI", "Father_BMI" = "Paternal BMI", "BMIz_1yr" = "Infant BMIz (1 year)", "Site_Vancouver" = "Study Site (Vancouver vs. Other)", "BF_3m_status_new" = "Exclusive Breastfeeding (3 months)", "BF_duration_imp" = "Breastfeeding Duration", "Shannon_bacteria_3m" = "Bacterial Alpha Diversity (3 months)", "Shannon_bacteria_12m" = "Bacterial Alpha Diversity (12 months)", "Prebirth_abs_oralIV_yn" = "Prenatal Antibiotics", "csec" = "Birth Mode", "Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)", "Mother_abs_birth_yn" = "Intrapartum Antibiotics", "sex" = "Infant Sex", "AS_bev4" = "Maternal AS Beverages in Pregnancy", "PCoA1_bacteria_3m" = "Bacterial Beta Diversity (3 months)", "PCoA1_bacteria_12m" = "Bacterial Beta Diversity (12 months)", "solids_6m" = "Introduction of Solids (6 months)"))
random_forest_fig

# Load packages for logistic regression
library(stats)
library(performance)

# Using dataframe from random forest, adjust for logistic regression by scaling numeric values to improve confidence interval visualization
y <- x
y$BF_duration_months <- y$BF_duration_imp/12
y$mom_bmi_best_5 <- y$mom_bmi_best/5
y$hei2010_5 <- y$hei2010/5
y$Father_BMI_5 <- y$Father_BMI/5

# Check dataframe prior to plotting
summary(y)

# Perform logistic regression to determine direction and significance of association between covariates and Chao1 pattern
mylogit <- glm(Chao1_dif_cat ~ sex + csec + Site_Vancouver + BF_3m_status_new + BF_duration_months + solids_6m + Mother_abs_birth_yn + Child_6moto1Y_abs + Prebirth_abs_oralIV_yn + mom_bmi_best_5 + hei2010_5 + AS_bev4 + Father_BMI_5 + BMIz_3mo + BMIz_1yr, data=y, family = "binomial")
summary(mylogit)

# Ensure there is no substantial collinearity in model
check_collinearity(mylogit)

# Determine confidence intervals for logistic regression results and create dataframe
confint(mylogit)
logisticr_df <- exp(cbind(OR = coef(mylogit), confint(mylogit)))
logisticr_df <- as.data.frame(logisticr_df)

# Added p-values, column headers & categories manually in Excel from logistic regression results output
#write.csv(logisticr_df, "logisticr_df.csv")

# Reload dataframe for to plot logisitic regression results
logisticr_df <- read.csv("logisticr_df.csv")

# Create column with stars to denote significant p-values
logisticr_df$stars <- cut(logisticr_df$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", "~",""))

# Order variables based on random forest mean decreasing Gini index to enable consistent plots
logisticr_df$Variable <- factor(logisticr_df$Variable, levels = c("Prebirth_abs_oralIV_ynYes", "BF_3m_status_newExclusive", "AS_bev43", "csecC-Section", "Child_6moto1Y_abs1", "sexFemale", "Site_VancouverOther", "Mother_abs_birth_ynYes", "solids_6mYes", "BF_duration_months", "BMIz_3mo", "hei2010_5", "mom_bmi_best_5", "Father_BMI_5", "BMIz_1yr"))

# Forest plot of logistic regression results (Figure S3)
log_regression_fp <- ggplot(logisticr_df, aes(x=Variable, y=OR, ymin=Lower, ymax=Upper, color = Factor))+
  geom_pointrange(position = position_dodge(width = 0.60))+
  geom_hline(yintercept = 1, linetype=2)+
  coord_flip()+
  xlab('')+
  scale_y_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)))+
  ylab("Adjusted OR (95% CI)")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper),width=0.1, position = position_dodge(width = 0.60))+
  geom_text(aes(x = Variable, y = OR, label=stars), color="black", size=6, hjust = 0.5, vjust = 0.1) + 
  theme_bw()+
  theme(plot.title = element_text(size=10, hjust = 0.5, face = "bold"),
        panel.background = element_blank(),
        axis.title=element_text(face = "bold", size = 10),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9),
        legend.text.align = 0,
        legend.position = "right")+
  scale_colour_manual(values=c("purple4","mediumorchid3","deepskyblue4","#35b779","olivedrab3"))+
  scale_x_discrete(labels=c("hei2010_5" = "Maternal Healthy Eating Index (per 5 units)", "BMIz_3mo" = "Infant BMIz (3 months)", "mom_bmi_best_5" = "Maternal BMI (per 5 units)", "Father_BMI_5" = "Paternal BMI (per 5 units)", "BMIz_1yr" = "Infant BMIz (1 year)", "Site_VancouverOther" = "Study Site (Vancouver vs. Other)", "BF_duration_months" = "Breastfeeding Duration (months)", "BF_3m_status_newExclusive" = "Exclusive Breastfeeding (3 months)", "solids_6mYes" = "Introduction of Solids (6 months)", "Prebirth_abs_oralIV_ynYes" = "Prenatal Antibiotics", "csecC-Section" = "Birth Mode (Vaginal vs. C-Section)", "Child_6moto1Y_abs1" = "Antibiotic Exposure (6-12 months)", "Mother_abs_birth_ynYes" = "Intrapartum Antibiotics", "sexFemale" = "Infant Sex (Male vs. Female)", "AS_bev43" = "Maternal AS Beverages in Pregnancy"))
log_regression_fp

