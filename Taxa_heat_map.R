#Taxa heat map for BMIz (Figure 3 and S6)

#load packages 
library(tidyverse)
library(phyloseq)
library(zCompositions)
library(CoDaSeq)
library(lmPerm)

#set working directory
setwd("/Users/mwgutierrez/Desktop/CHILD")

#load phyloseq object
load("ps_fungi")

#relativise and set the threshhold for taxa with greater than 2% mean relative abundance
d <- tax_glom(ps_fungi_clean,taxrank = "Genus") %>%
  transform_sample_counts(function(x) x*100 / sum(x)) %>%
  filter_taxa(function(x) mean(x) > 2, TRUE) 
d

#filter for 3 month samples
d_3M <- d %>%
  subset_samples(age == "three_months")

#generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(d_3M))
                             [c("Genus")]))  # to distinguish from "_" within tax ranks

#turn the otu_table into a data.frame
otu_export <- as.data.frame(otu_table(d_3M))
tmp <- names(otu_export)

#paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i])
}

#overwrite and clean up names
names(otu_export) <- names(tmp)
otu_export

otu_export2 <- otu_export %>% 
  rename_all(.funs = funs(sub("\\ ASV.*", "", names(otu_export)))) 

otu_export2 <- otu_export2%>%
  rename_all(.funs = funs(sub("\\g__", "", names(otu_export2)))) 
otu_export2

#extract dataframe 
sdf <- as(sample_data(d_3M), "data.frame")
sdf

#recode BF status to Zero/Partial vs Exclusive
sdf <- sdf %>%
  mutate(BF_3m_status=dplyr::recode(BF_3m_status,
                                    "Partial"="Zero/Partial",
                                    "Zero"="Zero/Partial",
                                    "Exclusive"="Exclusive"))
sdf$BF_3m_status <- as.factor(sdf$BF_3m_status)
sdf$BF_3m_status<- relevel(sdf$BF_3m_status, ref = "Exclusive")

##zero corrections and CLR transformations
#replace 0 values with an estimate of the probability that the zero is not 0
d.n0 <- cmultRepl(otu_export2,  label=0, method="CZM") 
#CLR transformation
d.n0.clr <- as.data.frame(codaSeq.clr(d.n0, samples.by.row=TRUE)) 

#check that otu table and metadata are in the exact same order (and have the same sampleIDs as rownames)
identical(rownames(d.n0.clr), rownames(sdf))

##model the effect of fungi at 3M (with covariates) on BMIz at 3M
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_3mo ~ x + sdf$BF_3m_status + sdf$hei2010 + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func

#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 3M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p) %>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_3M_3M <- m_coeff %>%
  filter(Factor == "Fungi 3M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_3M_3M

##model the effect of fungi at 3M (with covariates) on BMIz at 1y
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_1yr ~ x + sdf$BF_3m_status + sdf$hei2010 + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func

#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 3M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p)%>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_3M_1y <- m_coeff %>%
  filter(Factor == "Fungi 3M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_3M_1y

##model the effect of fungi at 3M (with covariates) on BMIz at 3y
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_3yr ~ x + sdf$hei2010 + sdf$BF_3m_status + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func

#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 3M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p)%>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_3M_3y <- m_coeff %>%
  filter(Factor == "Fungi 3M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_3M_3y

##model the effect of fungi at 3M (with covariates) on BMIz at 5y
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_5yr ~ x +  sdf$BF_3m_status + sdf$hei2010 + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func

#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 3M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p) %>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_3M_5y <- m_coeff %>%
  filter(Factor == "Fungi 3M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_3M_5y

#filter for 12 month samples
d_12M <- d %>%
  subset_samples(age == "twelve_months")

#generate a vector containing the full taxonomy path for all OTUs
wholetax <- do.call(paste, c(as.data.frame(tax_table(d_12M))
                             [c("Genus")]))  # to distinguish from "_" within tax ranks

#turn the otu_table into a data.frame
otu_export <- as.data.frame(otu_table(d_12M))
tmp <- names(otu_export)

#paste wholetax and OTU_ids together
for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i])
}

#overwrite and clean up names
names(otu_export) <- names(tmp)
otu_export

otu_export2 <- otu_export %>% 
  rename_all(.funs = funs(sub("\\ ASV.*", "", names(otu_export)))) 

otu_export2 <- otu_export2%>%
  rename_all(.funs = funs(sub("\\g__", "", names(otu_export2)))) 
otu_export2

#extract dataframe from phyloseq object
sdf <- as(sample_data(d_12M), "data.frame")
sdf

#recode BF status to Zero/Partial vs Exclusive
sdf <- sdf %>%
  mutate(BF_3m_status=dplyr::recode(BF_3m_status,
                                    "Partial"="Zero/Partial",
                                    "Zero"="Zero/Partial",
                                    "Exclusive"="Exclusive"))
sdf$BF_3m_status <- as.factor(sdf$BF_3m_status)
sdf$BF_3m_status<- relevel(sdf$BF_3m_status, ref = "Exclusive")
sdf

#zero corrections and CLR transformations
#replace 0 values with an estimate of the probability that the zero is not 0
d.n0 <- cmultRepl(otu_export2,  label=0, method="CZM") 

#CLR transformation
d.n0.clr <- as.data.frame(codaSeq.clr(d.n0, samples.by.row=TRUE)) 
 
#check that otu table and metadata are in the exact same order (and have the same sampleIDs as rownames)
identical(rownames(d.n0.clr), rownames(sdf))

#model the effect of fungi at 12M (with covariates) on BMI 1y 
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_1yr ~ x + sdf$BF_3m_status + sdf$hei2010 + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func

#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 12M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p)%>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_12M_1y <- m_coeff %>%
  filter(Factor == "Fungi 12M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_12M_1y

#model the effect of fungi at 12M (with covariates) on BMIz 3yr
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_3yr ~ x + sdf$BF_3m_status + sdf$hei2010 + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func
 
#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 12M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p)%>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_12M_3y <- m_coeff %>%
  filter(Factor == "Fungi 12M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_12M_3y

##model the effect of fungi at 12M (with covariates) on BMIz 5yrs
#define the function
lmp_func <- function(x) {lmp(sdf$BMIz_5yr ~ x + sdf$BF_3m_status + sdf$hei2010 + sdf$mom_bmi_best + sdf$Child_6moto1Y_abs, maxIter=10000, Ca = 0)} 
lmp_func

#use the lmp function on each taxa individually
set.seed(9999)
lmp_test <- lapply(d.n0.clr, FUN = lmp_func)

##get the summary statistics
#create summary statistics function
lmp_summ <- lapply(lmp_test, FUN= function(x) {summary(x)}) 

#CI intervals
lmp_confits <- as.data.frame(lapply(lmp_test, FUN= function(x) {confint(x)}))
lmp_confits

#get summary of coefficients
lmp_coeff <- as.data.frame(lapply(lmp_summ, function(x) {x["coefficients"]}))
lmp_coeff

m_coeff <- lmp_coeff %>%
  rownames_to_column("Factor")%>%
  filter(Factor != "(Intercept)")%>%
  mutate(Factor=dplyr::recode(Factor,
                              "x" = "Fungi 12M",
                              "sdf$BF_3m_status"="Zero/Partial vs. Exclusive BF",
                              "sdf$hei2010" = "Maternal HEI",
                              "sdf$mom_bmi_best" = "Maternal BMI",
                              "sdf$Child_6moto1Y_abs" = "Antibiotic Exposure (6-12 months)")) %>%
  gather(colnames(lmp_coeff), key = variable, value = value)%>%
  filter(!grepl(".Iter",variable))%>%
  separate(variable, c("ASV", "sequence", "extra", "variable", "p"), sep = "\\.")%>%
  dplyr::select(-extra,-p)%>%
  spread(variable, value) %>%
  group_by(Factor, ASV)%>%
  summarise_all(funs(trimws(paste(., collapse = ''))))
m_coeff

#adjust p-value and save results as a dataframe
m_coeff$Pr <- as.numeric(m_coeff$Pr)
m_coeff$Estimate <- as.numeric(m_coeff$Estimate)

results_12M_5y <- m_coeff %>%
  filter(Factor == "Fungi 12M") %>%
  mutate(p.adjust=p.adjust(Pr, method = "BH")) 
results_12M_5y

##make a large heat map with every BMIz timepoint 
#add new column to data frames for BMIz timepoint 
results_3M_3M_Time <- cbind(results_3M_3M, time = "3 months") 
results_3M_1y_Time <- cbind(results_3M_1y, time = "1 year") 
results_3M_3y_Time <- cbind(results_3M_3y, time = "3 years") 
results_3M_5y_Time <- cbind(results_3M_5y, time = "5 years") 
results_12M_1y_Time <- cbind(results_12M_1y, time = "1 year") 
results_12M_3y_Time <- cbind(results_12M_3y, time = "3 years") 
results_12M_5y_Time <- cbind(results_12M_5y, time = "5 years") 

results_all_time <- rbind(results_3M_3M_Time, results_3M_1y_Time, results_3M_3y_Time, results_3M_5y_Time, results_12M_1y_Time, results_12M_3y_Time, results_12M_5y_Time)


#create variable of stars based on p-values
results_all_time$stars <- cut(results_all_time$p.adjust, breaks=c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", "~", ""))  # Create column of significance labels

#get the min and max of the Estimate
m <- results_all_time %>%
  summarize(min=min(Estimate), max=max(Estimate))
m <- m[-1,]#if having issues with it grouping the min and max by Factor

Heatmap_limits <- c(m$min, m$max)

#plot heat map and save
results_all_time$time <- factor(results_all_time$time,levels = c("3 months", "1 year", "3 years", "5 years"))
results_all_time$Factor <- factor(results_all_time$Factor,levels = c("Fungi 3M", "Fungi 12M"))

jpeg(file = "/Users/mwgutierrez/Desktop/heatmap_genus.jpeg", width = 1300, height = 1500, units = "px", res = 300)   
Fig <- ggplot(results_all_time, aes(Factor, ASV, fill=Estimate)) +
  facet_wrap(~time, nrow=1, scales = "free_x") +
  geom_tile(height=0.8, width=0.8) +
  scale_fill_gradient2(low="#482878", mid="white", high="#35b779", limits=Heatmap_limits, breaks = c(-0.10, -0.05, 0, 0.05, 0.10, 0.15)) +
  theme_bw()+
  labs(x = "", y = "", fill="Beta coefficient") + 
  theme(axis.ticks = element_blank(), 
        axis.text.x=element_text(size=8, angle=45, vjust=1, hjust=1, color = "black"),
        axis.text.y=element_text(size=8, 
                                 color = "black", face = "italic", angle = 0),
        panel.grid.major=element_blank(), 
        legend.position = "bottom",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        strip.text.x = element_text(size = 8, face = "bold"), 
        legend.title=element_text(size=8, face = "bold"), 
        legend.key.width=unit(1.15,"cm"),
        legend.title.align = 0.5) +
  geom_text(aes(label=stars), color="black", size=5) +
  scale_y_discrete(limits = rev)
Fig
dev.off()
