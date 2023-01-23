library(phyloseq)
library(tidyverse)
library(lavaan)
library(vegan)
library(MVN); packageVersion("MVN")

theme_set(theme_bw())

load('ps_fungi_clean_CHILD_21_May_2021')
ps_fungi_clean

# Porting Data to vegan function
veganotu = function(data) {
  require("vegan")
  OTU = otu_table(data)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

sdf <- as(sample_data(ps_fungi_clean), "data.frame")

tax <- as.data.frame(as(tax_table(ps_fungi_clean), "matrix"))

sort(unique(tax$Genus))

select_genus_3m <- c( "g__Candida","g__Cladosporium", "g__Fomitopsis",   "g__Malassezia",  "g__Mycosphaerella", "g__Rhodotorula" )

genus_unique_in_3 <- c( "g__Cladosporium", "g__Fomitopsis" )

select_genus_12m <- c( "g__Candida", "g__Malassezia", "g__Mycosphaerella", "g__Rhodotorula",  "g__Saccharomyces")

genus_unique_in_12 <- c( "g__Saccharomyces" )


d_genus <- tax_glom(ps_fungi_clean,taxrank = "Genus") %>%
  transform_sample_counts(function(x) x*100/sum(x)) %>% 
  subset_taxa(Genus == "g__Candida")
d_genus

d_fung <- as.data.frame(veganotu(d_genus))

d_fung <- as.data.frame(veganotu(d_genus))%>%
  rename(g__Candida=TGTGGTGGTGGGTCCTCCGCTTATTGATATGCTTAAGTTCAGCGGGTAGTCCTACCTGATTTGAGGTCGAATTTGGAAGAAGTTTTGGAGTTTGTACCAATGAGTGGAAAAAACCTATCCATTAGTTTATACTCCGCCTTTCTTTCAAGCAAACCCAGCGTATCGCTCAACACCAAACCCGAGGGTTTGAGGGAGAAATGACGCTCAAACAGGCATGCCCTTTGGAATACCAAAGGGCGCAATGTGCGTTCAAAGATTCGATGACTCACGGCGGCCTGACG)
head(d_fung)

identical(rownames(sdf) , rownames(d_fung))

sdf_genus <- cbind(sdf, d_fung)%>%
  select(CHILD_ID1, g__Candida, sex, BMIz_3mo , BMIz_1y.1,  BMIz_3y.1,  BMIz_5y, Chao1_bacteria, PCoA1_bacteria, Chao1_dif_cat, PCoA1_fungi_novst, mom_bmi_best, hei2010, BF_3m_status, Father_BMI, Sample_time, Chao1)%>%
  mutate(BM=recode(BF_3m_status, 
                   "Zero"="0",
                   "Partial"="0",
                   "Exclusive"="1"))%>%
  mutate(Sex=recode(sex, 
                    "Male"="0",
                    "Female"="1"))

sdf_genus$Sex <- as.numeric(sdf_genus$Sex)

head(sdf_genus)

d3 <- sdf_genus%>% 
  filter(Sample_time == "3")%>%
  arrange(CHILD_ID1)%>%
  select(CHILD_ID1, BMIz_5y, BMIz_3y.1, BMIz_1y.1, BMIz_3mo, g__Candida, PCoA1_bacteria, mom_bmi_best, BM, hei2010, Sex , Chao1, Chao1_dif_cat, PCoA1_fungi_novst)%>%
  rename(Candida_3m = g__Candida, Bacteria_3m = PCoA1_bacteria, Chao1_3m = Chao1, fungi_pco1_3m = PCoA1_fungi_novst)

head(d3)

d12 <- sdf_genus%>% 
  filter(Sample_time == "12")%>%
  arrange(CHILD_ID1) %>%
  select(CHILD_ID1, g__Candida, PCoA1_bacteria, Chao1, PCoA1_fungi_novst)%>%
  rename(Candida_1y = g__Candida, Bacteria_1y = PCoA1_bacteria, Chao1_1y = Chao1, fungi_pco1_1y = PCoA1_fungi_novst)
head(d12)

d_all <- merge(d3, d12, by="CHILD_ID1")%>%
  filter(complete.cases(.)) %>%
  column_to_rownames("CHILD_ID1")
head(d_all)

d_all$BM <- as.numeric(d_all$BM)
d_all$Chao1_3m  <- as.numeric(d_all$Chao1_3m)
d_all$Chao1_1y <- as.numeric(d_all$Chao1_1y)


# Candida at 3 and 12 months

d_model_1 <- d_all %>%
  mutate(Candida_3m = log10(Candida_3m), Candida_1y=log10(Candida_1y))%>%
  select(-Chao1_1y, -  Chao1_3m ,  -Chao1_dif_cat, -fungi_pco1_3m, -fungi_pco1_1y)
head(d_model_1)

set.seed(5879)
result <- mvn(data = d_model_1,  mvnTest = "mardia", univariateTest = "AD", univariatePlot = "histogram",
              multivariatePlot = "qq")

result$multivariateNormality

set.seed(1894068)

model1 <- ' # regressions
            
              BMIz_5y ~ BMIz_3y.1   + Candida_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3y.1  ~ BMIz_1y.1  + Candida_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_1y.1 ~ BMIz_3mo + Candida_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3mo ~ Candida_3m  + Bacteria_3m + mom_bmi_best + BM + hei2010 + Sex 
              
              # define the association of bacteria and fungi
              Candida_1y    ~ Candida_3m + mom_bmi_best + BM 
              Bacteria_1y  ~  mom_bmi_best + BM + Bacteria_3m 
              
             
     #Chao1_bacteria ~~  0* mom_bmi_best  
     
     #BM ~~  site  
     mom_bmi_best ~~ 0*BM  
     hei2010 ~~ 0*Sex  
     #g__Candida    ~~ Chao1_bacteria 
     #BMIz_5y ~~ 0* Site   
     Candida_3m ~~  Bacteria_3m
     Candida_1y ~~  Bacteria_1y
      
fit1 <- sem(model1, data=d_model_1, std.lv = TRUE, orthogonal = TRUE, std.ov = TRUE, fixed.x = FALSE, se = "bootstrap", verbose = FALSE, bootstrap = 1000)
summary(fit1, fit.measures=TRUE, standardized=TRUE, rsq = TRUE)

#tiff(file = "figures/38_sem1.tiff", width = 2200, height = 900, units = "px", res = 300)
#semPaths(fit1, "std",  edge.label.cex = 0.85, curvePivot = FALSE, exoVar = FALSE, fade = FALSE, curve = TRUE)
#dev.off()

fitMeasures(fit1, c("cfi", "rmsea", "srmr"))
AIC(fit1)
BIC(fit1)

# Fungi Chao1 (Richness) at 3 and 12 months

d_model_2 <- d_all %>%
  #mutate(Candida_3m = log10(Candida_3m), Candida_1y=log10(Candida_1y))%>%
  select(-Candida_3m, -Candida_1y,  -Chao1_dif_cat, -fungi_pco1_3m, -fungi_pco1_1y)
head(d_model_2)

set.seed(5879)
result <- mvn(data = d_model_2,  mvnTest = "mardia", univariateTest = "AD", univariatePlot = "histogram",
              multivariatePlot = "qq")

result$multivariateNormality

set.seed(1894068)

model1 <- ' # regressions
            
              
              BMIz_5y ~ BMIz_3y.1  + Chao1_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3y.1  ~ BMIz_1y.1 +  Chao1_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_1y.1 ~ BMIz_3mo + Chao1_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3mo ~ Chao1_3m  + Bacteria_3m + mom_bmi_best + BM + hei2010 + Sex 
              
              
              Chao1_1y    ~ Chao1_3m + mom_bmi_best + BM 
              Bacteria_1y  ~  mom_bmi_best + BM + Bacteria_3m   
             
     #Chao1_bacteria ~~  0* mom_bmi_best  
     
     #BM ~~  site  
     mom_bmi_best ~~ 0*BM  
     hei2010 ~~ 0*Sex  
     #g__Candida    ~~ Chao1_bacteria 
     #BMIz_5y ~~ 0* Site   
     Chao1_3m ~~  Bacteria_3m
     Chao1_1y ~~  Bacteria_1y
      
               '
fit1 <- sem(model1, data=d_model_2, std.lv = TRUE, orthogonal = TRUE, std.ov = TRUE, fixed.x = FALSE, se = "bootstrap", verbose = FALSE, bootstrap = 1000)
summary(fit1, fit.measures=TRUE, standardized=TRUE, rsq = TRUE)

fitMeasures(fit1, c("cfi", "rmsea", "srmr"))
AIC(fit1)
BIC(fit1)

# Fungi PcoA1 (Beta Diversity) at 3 and 12 months

d_model_3 <- d_all %>%
  select(-Candida_3m, -Candida_1y,  -Chao1_dif_cat, -Chao1_3m, -Chao1_1y)
head(d_model_3)

set.seed(5879)
result <- mvn(data = d_model_3,  mvnTest = "mardia", univariateTest = "AD", univariatePlot = "histogram",
              multivariatePlot = "qq")

result$multivariateNormality

set.seed(1894068)

model1 <- ' # regressions
            
              
              BMIz_5y ~ BMIz_3y.1  + fungi_pco1_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3y.1  ~ BMIz_1y.1 +  fungi_pco1_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_1y.1 ~ BMIz_3mo + fungi_pco1_1y  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3mo ~ fungi_pco1_3m  + Bacteria_3m + mom_bmi_best + BM + hei2010 + Sex 
              
              
              fungi_pco1_1y    ~ fungi_pco1_3m + mom_bmi_best + BM  
              Bacteria_1y  ~  mom_bmi_best + BM + Bacteria_3m   
             
     #Chao1_bacteria ~~  0* mom_bmi_best  
     
     #BM ~~  site  
     mom_bmi_best ~~ 0*BM  
     hei2010 ~~ 0*Sex  
     #g__Candida    ~~ Chao1_bacteria 
     #BMIz_5y ~~ 0* Site   
     fungi_pco1_3m ~~  Bacteria_3m
     fungi_pco1_1y ~~  Bacteria_1y
      
               '
fit1 <- sem(model1, data=d_model_3, std.lv = TRUE, orthogonal = TRUE, std.ov = TRUE, fixed.x = FALSE, se = "bootstrap", verbose = FALSE, bootstrap = 1000)
summary(fit1, fit.measures=TRUE, standardized=TRUE, rsq = TRUE)

fitMeasures(fit1, c("cfi", "rmsea", "srmr"))
AIC(fit1)
BIC(fit1)

# Fungi Richness (Chao1) Pattern

m <- d_all %>% filter(Chao1_dif_cat != "Unchanged")%>%
  mutate(Chao1_dif_cat=recode(Chao1_dif_cat, 
                              "Decrease"="0",
                              "Increase"="1"))
m

m$Chao1_dif_cat <- as.numeric(m$Chao1_dif_cat)

d_model_4 <- m %>%
  select(-Candida_3m, -Candida_1y,  -Chao1_3m, -Chao1_1y, -fungi_pco1_3m, -fungi_pco1_1y)
head(d_model_4)

set.seed(5879)
result <- mvn(data = d_model_4,  mvnTest = "mardia", univariateTest = "AD", univariatePlot = "histogram",
              multivariatePlot = "qq")

result$multivariateNormality

set.seed(1894068)

model1 <- ' # regressions
            
              
              BMIz_5y ~ BMIz_3y.1  + Chao1_dif_cat  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3y.1  ~ BMIz_1y.1 +  Chao1_dif_cat  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_1y.1 ~ BMIz_3mo + Chao1_dif_cat  + Bacteria_1y + mom_bmi_best + BM + hei2010 + Sex 
              
              BMIz_3mo ~ Chao1_dif_cat  + Bacteria_3m + mom_bmi_best + BM + hei2010 + Sex 
              
              
              Chao1_dif_cat    ~ Bacteria_3m + Bacteria_1y + mom_bmi_best + BM  
              Bacteria_1y  ~  mom_bmi_best + BM + Bacteria_3m 
                     
     #Chao1_bacteria ~~  0* mom_bmi_best  
     
     #BM ~~  site  
     mom_bmi_best ~~ 0*BM  
     hei2010 ~~ 0*Sex  
     #g__Candida    ~~ Chao1_bacteria 
     #BMIz_5y ~~ 0* Site   
     #Chao1_3m ~~  Bacteria_3m
     #Chao1_1y ~~  Bacteria_1y
      
 
               '
fit1 <- sem(model1, data=d_model_4, std.lv = TRUE, orthogonal = TRUE, std.ov = TRUE, fixed.x = FALSE, se = "bootstrap", verbose = FALSE, bootstrap = 1000)
summary(fit1, fit.measures=TRUE, standardized=TRUE, rsq = TRUE)

fitMeasures(fit1, c("cfi", "rmsea", "srmr"))
AIC(fit1)
BIC(fit1)
