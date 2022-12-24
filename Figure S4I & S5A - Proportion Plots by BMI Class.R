#Proportion plots for BMI classes (Figure S4 and S5)

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

#remove the participants with unchanged or NA Chao1 pattern for this analysis
d_filt <- d_filt %>%
  filter(Chao1_dif_cat != "NA") %>%
  filter(Chao1_dif_cat != "Unchanged")
d_filt

#create long format based on BMI class variables
d_filt_long_class<- d_filt%>%
  gather(BMI_birth_class, BMI_3M_class,	BMI_1y_class,	BMI_3y_class, BMI_5y_class, key = variable, value = BMI_class)%>%
  separate(variable, into = c("first", "third"), sep = -9)%>%
  dplyr::select(-first)%>%
  mutate(Time=dplyr::recode(third,
                            "rth_class" = "Birth",
                            "_3M_class"="3 months",
                            "_1y_class"="1 year",
                            "_3y_class"="3 years",
                            "_5y_class" = "5 years"))%>%
  dplyr::select(-third)
d_filt_long_class

#prepare data for plotting
d_filt_long_class$Time <- factor(d_filt_long_class$Time,levels = c("Birth", "3 months", "1 year", "3 years", "5 years"))

d_filt_long_class$BMI_class <- factor(d_filt_long_class$BMI_class,levels = c("Risk of overweight/Overweight", "Normal"))

d_filt_long_class_noNA <- d_filt_long_class %>%
  filter(BMI_class != "NA")

prop_sdf <- d_filt_long_class_noNA %>%
  group_by(Time, Chao1_dif_cat, BMI_class)%>%
  dplyr::summarise(n = n()) %>%
  mutate(freq = n*100 / sum(n))
prop_sdf

#plot proportion of BMI classes based on Chao1 pattern stratified by time point and save
jpeg(file = "prop_BMIz_class_age.jpeg", width = 2200, height = 600, units = "px", res = 300)
fig <- ggplot(prop_sdf, aes(x = Chao1_dif_cat, y = freq, fill = BMI_class)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  facet_wrap(~Time, scales = "free_x", nrow = 1) +
  labs(y = "% BMIz class", x = "Richness Pattern", fill = "BMIz class") +
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(labels = c("Risk of Overweight/ \nOverweight", "Normal"), values = c("#3e4989", "#26828e"))+
  theme(legend.position="right") + 
  theme(panel.spacing.y = unit(0.6, "lines"))+
  theme(panel.background = element_blank(),
        axis.title = element_text(face = "bold", size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 8), 
        legend.text = element_text(size = 7),
        legend.text.align = 0,
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 8))

fig
dev.off()

#perform chi-squared test (repeat for each time point)
t <- d_filt_long_class_noNA %>%
  filter(Time == "5 years") 

tbl = table(t$BMI_class, t$Chao1_dif_cat)

chisq.test(tbl)
