#Regression interactions (Figure 2, Table S3 and Figure S5)

#load packages 
library(tidyverse)
library(stats)
library(stargazer)
library(effects)

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

#create long format based on BMI variables
d_filt_long <- d_filt %>%
  gather(BMIz_birth, BMIz_3mo,	BMIz_1yr,	BMIz_3yr, BMIz_5yr, key = variable, value = BMI_z) %>%
  separate(variable, into = c("first", "third"), sep = -4) %>%
  dplyr::select(-first) %>%
  mutate(Time=dplyr::recode(third, 
                            "irth" = "0",
                            "_3mo"="90",
                            "_1yr"="365",
                            "_3yr"="1095",
                            "_5yr" = "1825")) %>%
  dplyr::select(-third)
d_filt_long$Time <- as.numeric(d_filt_long$Time)
d_filt_long


#recode breastfeeding status to Zero/Partial vs Exclusive and sex variable
d_filt_long <- d_filt_long %>%
  mutate(BF_3m_status=dplyr::recode(BF_3m_status,
                                    "Partial"="Zero/Partial",
                                    "Zero"="Zero/Partial",
                                    "Exclusive"="Exclusive")) %>%
  mutate(sex=dplyr::recode(sex,
                           "1"="Male",
                           "0"="Female"))
d_filt_long

#create models for BMIz outcomes where Chao1 pattern is independent or interacting with the variable of interest
mod1 <- lm(BMI_z ~ poly(Time,2) + mom_bmi_best+Chao1_dif_cat , data=d_filt_long)
mod2 <- lm(BMI_z ~ poly(Time,2) + mom_bmi_best*Chao1_dif_cat , data=d_filt_long)

#determine the main effects and interactions
stargazer(mod1, mod2, type="text", 
          column.labels = c("Main Effects", "Interaction"), 
          intercept.bottom = FALSE, 
          single.row=FALSE,     
          notes.append = FALSE, 
          header=FALSE,
          report=('vc*p')) 

#prepare data for plotting
m <- effect('mom_bmi_best*Chao1_dif_cat', mod2,
            xlevels=list(Chao1_dif_cat = c("Increase","Decrease"),
                         mom_bmi_best = c(18, 25, 30, 35)),
            se=TRUE, confidence.level=.95, typical=mean)

#put data in data frame 
m <- as.data.frame(m)
head(m)

#create a factor of the Chao1 pattern variable used in the interaction                   
m$Chao1_dif_cat <- factor(m$Chao1_dif_cat,
                          levels=c("Decrease","Increase"))

#create a factor of the Maternal BMI variable used in the interaction 
m$mom_bmi_best <- factor(m$mom_bmi_best,
                         levels=c(18, 25, 30, 35),
                         labels=c("18-25", "25-30", "30-35", ">35"))

#plot and save 
jpeg(file = "interaction_momBMI.jpeg", width = 950, height = 560, units = "px", res = 300)       
fig<-ggplot(data=m, aes(x=mom_bmi_best, y=fit, group=Chao1_dif_cat))+
  geom_line(size=1, aes(color=Chao1_dif_cat))+
  ylim(0,2)+
  ylab("BMIz (birth-5yrs)")+
  xlab("Maternal BMI")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=5))+
  scale_color_manual(name = "Richness pattern", values=c("#2db27d","#481b6d"))+
  theme_bw()+
  theme(plot.title = element_text(size=8, hjust = 0.5, face = "bold"))+
  theme(axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 8), 
        legend.text = element_text(size = 7),
        legend.text.align = 0)
fig
dev.off()
