#Boxplots for study variables based on Chao1 pattern (Figure S4)

#load packages 
library(tidyverse)

#set working directory
setwd("/Users/mackenziegutierrez/Desktop/CHILD")

#load phyloseq object
load("ps_fungi")

#extract sample data from phyloseq object to use for analysis
df <- as(sample_data(ps_fungi_clean), "data.frame")

#filter based on sample time to remove duplicate data 
d_filt <- df %>%
  filter(Sample_time == "3")

#Remove the participants with unchanged or NA Chao1 pattern for this analysis
d_filt <- d_filt %>%
  filter(Chao1_dif_cat != "NA") %>%
  filter(Chao1_dif_cat != "Unchanged")
d_filt

#plot variable of interest by Choa1 pattern and save
d_filt$Chao1_dif_cat <- factor(d_filt$Chao1_dif_cat,
                      levels=c("Decrease","Increase"))

jpeg(file = "boxplot_BMIzbirth.jpeg", width = 950, height = 560, units = "px", res = 300) 
Fig <- ggplot(d_filt, aes(x= Chao1_dif_cat, y = BMIz_birth, fill = Chao1_dif_cat )) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  geom_jitter(color="black", size=0.5, alpha=1) +
  scale_fill_manual(name = "Richness trend", values=c("#2db27d","#481b6d")) +
  labs(y = "BMIz", x = "Richness trend") +
  theme(axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 8), 
        legend.text = element_text(size = 7),
        legend.text.align = 0,
        plot.title = element_text(hjust = 0.5, face = "bold", size = 9)) +
  ggtitle("Birth") +
  scale_y_continuous(limits = c(-2.0, 4.0), labels=c(-2.0, 0.0, 2.0, 4.0))
Fig
dev.off()

#test for normality and homogeneity of variance
shapiro.test(d_filt$BMIz_birth)

var.test(BMIz_birth ~ Chao1_dif_cat, d_filt, 
         alternative = "two.sided")

#perform t-test
Chao1_increase <- subset(d_filt, d_filt$Chao1_dif_cat == "Increase")
Chao1_decrease <- subset(d_filt, d_filt$Chao1_dif_cat == "Decrease")
t.test(Chao1_increase$BMIz_birth, Chao1_decrease$BMIz_birth, paired=FALSE, var.equal = TRUE)