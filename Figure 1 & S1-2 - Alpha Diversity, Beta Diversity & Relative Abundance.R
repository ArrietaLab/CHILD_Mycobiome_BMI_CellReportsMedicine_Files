## Fungal Alpha Diversity, Beta Diversity & Relative Abundance (Figures 1 & S1-2)

# Load packages
library(tidyverse)
library(gpplot2)
library(phyloseq)

# Load phyloseq object
load("ps_fungi_clean")

# Extract sample data from phyloseq object to use for analysis 
fungi_data <- as(sample_data(ps_fungi_clean), "data.frame")

# Set seed for reproducible results
set.seed(60145)

## Alpha Diversity Overall (Figure 1A)

# Estimate fungal Shannon and Chao1 alpha diversity measures
adiv_fungi <- estimate_richness(ps_fungi_clean, measures = c("Shannon","Chao1"))

# Create new metadata columns for Chao1 and Shannon
fungi_data$Shannon <- adiv_fungi$Shannon
fungi_data$Chao1 <- adiv_fungi$Chao1

# Create filtered dataframes by timepoint
fungi_data_3 <- fungi_data %>%
  filter(Sample_time == "3")

fungi_data_12 <- fungi_data %>%
  filter(Sample_time == "12")

# Determine mean and standard deviation of fungal alpha diversity metrics at 3 and 12 months
mean(fungi_data_3$Shannon)
sd(fungi_data_3$Shannon)
mean(fungi_data_12$Shannon)
sd(fungi_data_12$Shannon)

mean(fungi_data_3$Chao1)
sd(fungi_data_3$Chao1)
mean(fungi_data_12$Chao1)
sd(fungi_data_12$Chao1)

# Plot Chao1 by sample time (Figure 1A) 
chao1_fig <- ggplot(fungi_data, aes(x = Sample_time, y = Chao1, fill = Sample_time)) + 
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 1) +
  geom_jitter(color="black", size=0.8, alpha=1) + 
  labs(x="Age (Months)",y="Richness (Chao1)")+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#228c8d","#443b84"))+
  theme_bw()+
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9),
        legend.position = "none")
chao1_fig

# Test for normality & perform Wilcox test
shapiro.test(fungi_data$Chao1)
wilcox.test(Chao1 ~ Sample_time, data = fungi_data)

# Plot Shannon by sample time (Figure 1A) 
shannon_fig <- ggplot(fungi_data, aes(x = Sample_time, y = Shannon, fill = Sample_time, color = Sample_time)) + 
  geom_boxplot(alpha = 0.6) +
  geom_jitter(aes(color = Sample_time), position = position_jitter(0.2),  size = 1.2) + 
  stat_summary(fun.data = summary_stats, size = 0.25)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(colour = NULL)+ 
  labs(x="Age (Months)")+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#228c8d","#443b84"))+
  theme_bw()+
  ylab("Alpha Diversity (Shannon)")+
  theme(strip.text.y = element_text(face = "italic"))+
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 9))+
  theme(legend.position = "none")
shannon_fig

# Test for normality & perform Wilcox test
shapiro.test(fungi_data$Shannon)
wilcox.test(Chao1 ~ Sample_time, data = fungi_data)

## Beta Diversity Overall (Figure 1B)

# Load packages for beta-diversity
library(DESeq2)
library(vegan)

# Create function geo means for variance stabilizing transformation
gm_mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}

# Variance stabilizing transformation 
vst_phyloseq_dseq_fungi <- phyloseq_to_deseq2(ps_fungi_clean, ~Sample_time)
vst_phyloseq_dseq_fungi = estimateSizeFactors(vst_phyloseq_dseq_fungi, geoMeans = apply(counts(vst_phyloseq_dseq_fungi), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(vst_phyloseq_dseq_fungi, blind = TRUE)
vst_blind_mat_fungi <- SummarizedExperiment::assay(vst_blind)
vst_blind_mat_fungi <- t(vst_blind_mat_fungi) 
vst_blind_mat_fungi[which(vst_blind_mat_fungi < 0)] <- 0 
dists <- dist(t(assay(vst_phyloseq_dseq_fungi)))

# Computing Bray-Curtis Dissimilarities and PCoA scores
comm_vst_blind_mat_fungi<- vegdist(vst_blind_mat_fungi, "bray")
PCoA_comm_vst_blind_mat_fungi <- capscale(comm_vst_blind_mat_fungi ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat_fungi$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat_fungi$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat_fungi)$sites

# Save PCoA scores into metadata tables
row.names(fungi_data) == row.names(scores(PCoA_comm_vst_blind_mat_fungi)$sites)
fungi_data$PCoA1_fungi <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,1]
fungi_data$PCoA2_fungi <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,2]

# PCoA plot by sample time (Figure 1B)
PCoA_fungi <- qplot(PCoA1_fungi, PCoA2_fungi, xlab = "PCoA1", ylab = "PCoA2",
                    size = 0.5, shape = Sample_time, fill = Sample_time, color = Sample_time, data = (fungi_data))

PCoA_fungi <- PCoA_fungi +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Sample_time, color = Sample_time)) +
  scale_shape_manual(name = "Age (Months)", values = c(21:24)) + 
  theme_bw() + 
  xlab("PCoA1 (10.7%)")+
  ylab("PCoA2 (5.2%)")+
  theme(legend.title = element_text(colour = "black", size = 9.5, face = "bold"),
        legend.text = element_text(colour = "black", size = 9.5),
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 10.5, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Age (Months)", values = c("black","black")) + 
  scale_fill_manual(name = "Age (Months)", values = c("#228c8d","#443b84")) + 
  guides(size = FALSE) 
PCoA_fungi

# PERMANOVA by sample time
permanova_fungi <- adonis(comm_vst_blind_mat_fungi ~ Sample_time, data = fungi_data, permutations = 999)
permanova_fungi$aov.tab

## Alpha Diversity Patterns (Figures 1C & S1A)
# Excel: Fungal Chao1 and Shannon alpha diversity metrics assessed individually by timepoint to determine the direction of change from 3 months to 1 year.
# Chao1 at 1 year - Chao1 at 3 months = Chao1 delta 
# Chao1 delta classified into increase, decrease, or unchanged categories
# Process repeated for Shannon index
# All subsequent analyses done using Chao1 and Shannon categories for change from 3 months to 1 year

# Reload metadata with alpha diversity patterns
fungi_data <- read.csv("fungi_data_adiv_patterns.csv")

# Rereate filtered dataframes by timepoint
fungi_data_3 <- fungi_data %>%
  filter(Sample_time == "3")

fungi_data_12 <- fungi_data %>%
  filter(Sample_time == "12")

# Filter data and relabel column values for Chao1 trend plot.
fungi_data_chao1 <- fungi_data %>% filter(Chao1_dif_cat != "NA")
fungi_data_chao1$Chao1_dif_cat_num <- fungi_data_chao1$Chao1_dif_cat # Relabel to include count per category
fungi_data_chao1 <- fungi_data_chao1 %>% mutate(Chao1_dif_cat_num = dplyr::recode(Chao1_dif_cat_num,
                                                                                "Increase" = "Increase (n = 25)",
                                                                                "Decrease" = "Decrease (n = 63)",
                                                                                "Unchanged" = "Unchanged (n = 3)"))

# Create filtered dataframes for Chao1 trend and remove NA values
fungi_data_3_chao1 <- fungi_data_3 %>%
  filter(Chao1_dif_cat != "NA")

fungi_data_12_chao1 <- fungi_data_12 %>%
  filter(Chao1_dif_cat != "NA")

# Plot Chao1 trends over the first year of life per individual (Figure 1C).
fig_chao1_pattern <- ggplot(fungi_data_chao1, aes(x = Sample_time, y = Chao1, group = CHILD_ID1, color = Chao1_dif_cat_num)) +
  geom_point() +
  geom_line(show.legend = F) +
  xlab("Age (Months)")+
  ylab("Richness (Chao1)")+
  labs(color = "Richness Pattern")+ 
  scale_color_manual(values = c("#2db27d", "#481b6d","goldenrod2"))+
  theme_bw()+
  theme(strip.text.y = element_text(face = "italic"),
        plot.title = element_text(size=9, hjust = 0.5),
        panel.background = element_blank(),
        axis.title=element_text(face = "bold", size = 9),
        axis.text = element_text(size = 8.5),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9),
        legend.text.align = 0)
fig_chao1_pattern

# Paired t-test for Chao1 trend from 3 to 12 months per individual
t.test(fungi_data_3_chao1$Chao1, 
       fungi_data_12_chao1$Chao1, 
       paired=TRUE, 
       conf.level=0.95)

# Filter data and relabel column values for Shannon trend plot.
fungi_data_shannon <- fungi_data %>% filter(Shannon_dif_cat != "NA")
fungi_data_shannon$Shannon_dif_cat_num <- fungi_data_shannon$Shannon_dif_cat
fungi_data_shannon <- fungi_data_shannon %>% mutate(Shannon_dif_cat_num = dplyr::recode(Shannon_dif_cat_num,
                                                                                  "Increase" = "Increase (n = 18)",
                                                                                  "Decrease" = "Decrease (n = 73)"))
# Create filtered dataframes for Chao1 trend and remove NA values
fungi_data_3_shannon <- fungi_data_3 %>%
  filter(Shannon_dif_cat != "NA")

fungi_data_12_shannon <- fungi_data_12 %>%
  filter(Shannon_dif_cat != "NA")

# Plot Shannon trends over the first year of life per individual (Figure S1A)
fig_shannon_pattern <- ggplot(fungi_data_shannon, aes(x = Sample_time, y = Shannon, group = CHILD_ID1, color = Shannon_dif_cat_num)) +
  geom_point() +
  geom_line(show.legend = F) +
  xlab("Age (Months)")+
  ylab("Alpha Diversity (Shannon)")+
  labs(color = "Alpha Diversity Pattern")+ 
  scale_color_manual(values = c("#2db27d", "#481b6d","goldenrod2"))+
  theme_bw()+
  theme(strip.text.y = element_text(face = "italic"),
        plot.title = element_text(size=9, hjust = 0.5),
        panel.background = element_blank(),
        axis.title=element_text(face = "bold", size = 10),
        axis.text = element_text(size = 9.5),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9),
        legend.text.align = 0)
fig_shannon_pattern

# Paired t-test for Shannon trend from 3 to 12 months per individual
t.test(fungi_data_3_shannon$Shannon, 
       fungi_data_12_shannon$Shannon, 
       paired=TRUE, 
       conf.level=0.95)

## Beta Diversity Patterns (Figures 1D & S1B)
# Repeat Bray-Curtis & PCoA plot using Chao1 and Shannon patterns

# Filter ps object by sample time and to remove unchanged/NA values for Chao1
ps_fungi_clean_chao1 <- subset_samples(ps_fungi_clean, Chao1_dif_cat != "NA")
ps_fungi_clean_chao1 <- subset_samples(ps_fungi_clean_chao1, Chao1_dif_cat != "Unchanged")
ps_fungi_clean_chao1_3 <- subset_samples(ps_fungi_clean_chao1, Sample_time == 3)

# Variance Stabilizing Transformation
vst_phyloseq_dseq_fungi <- phyloseq_to_deseq2(ps_fungi_clean_chao1_3, ~Chao1_dif_cat)

# Convert counts to integer
vst_phyloseq_dseq_fungi = estimateSizeFactors(vst_phyloseq_dseq_fungi, geoMeans = apply(counts(vst_phyloseq_dseq_fungi), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(vst_phyloseq_dseq_fungi, blind = TRUE)
vst_blind_mat_fungi <- SummarizedExperiment::assay(vst_blind) # Extract transformed OTU table
vst_blind_mat_fungi <- t(vst_blind_mat_fungi) # Transpose data
vst_blind_mat_fungi[which(vst_blind_mat_fungi < 0)] <- 0 # Set counts less than 0 to zero
dists <- dist(t(assay(vst_phyloseq_dseq_fungi)))

# Computing Bray-Curtis Dissimilarities and PCoA
comm_vst_blind_mat_fungi<- vegdist(vst_blind_mat_fungi, "bray")
PCoA_comm_vst_blind_mat_fungi <- capscale(comm_vst_blind_mat_fungi ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat_fungi$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat_fungi$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat_fungi)$sites

# Filter dataframe to remove unchanged Chao1 values
fungi_data_3_chao1_clean <- fungi_data_3_chao1 %>% filter(Chao1_dif_cat != "Unchanged")

# Save scores into metadata tables
row.names(fungi_data_3_chao1_clean) == row.names(scores(PCoA_comm_vst_blind_mat_fungi)$sites)
fungi_data_3_chao1_clean$PCoA1_fungi_3 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,1]
fungi_data_3_chao1_clean$PCoA2_fungi_3 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,2]

# Relabel column values for Chao1 trend beta-diversity plot
fungi_data_3_chao1_clean$Chao1_dif_cat_num <- fungi_data_3_chao1_clean$Chao1_dif_cat
fungi_data_3_chao1_clean <- fungi_data_3_chao1_clean %>% mutate(Chao1_dif_cat_num = dplyr::recode(Chao1_dif_cat_num,
                                                                                      "Increase" = "Increase (n = 25)",
                                                                                      "Decrease" = "Decrease (n = 63)"))

# PCoA plot for beta diversity at 3 months by Chao1 trend (Figure 1D)
PCoA_fungi_chao1_3 <- qplot(PCoA1_fungi_3, PCoA2_fungi_3, xlab = "PCoA1", ylab = "PCoA2",
                    size = 0.5, shape = Sample_time, fill = Chao1_dif_cat, color = Chao1_dif_cat, data = (fungi_data_3_chao1_clean))

PCoA_fungi_chao1_3 + stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Chao1_dif_cat, color = Chao1_dif_cat)) +
  scale_shape_manual(name = "Age (Months)", values = c(21:24)) +
  theme_bw() + 
  xlab("PCoA1 (6.3%)")+
  ylab("PCoA2 (4.8%)")+
  ggtitle("3 Months")+
  #ylim(-2,2)+
  #xlim(-2,3)+
  theme(legend.title = element_text(colour = "black", size = 9, face = "bold"),
        legend.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Richness Pattern", values = c("black","black")) + 
  scale_fill_manual(name = "Richness Pattern", values = c("#2db27d", "#481b6d")) + 
  guides(size = FALSE)

# PERMANOVA by Chao1 trend at 3 months
permanova_fungi <- adonis(comm_vst_blind_mat_fungi ~ Chao1_dif_cat, data = fungi_data_3_chao1_clean, permutations = 999)
permanova_fungi$aov.tab

# Repeat for Chao1 trend at 12 months

# Filter ps object by sample time
ps_fungi_clean_chao1_12 <- subset_samples(ps_fungi_clean_chao1, Sample_time == 12)

# Variance Stabilizing Transformation
vst_phyloseq_dseq_fungi <- phyloseq_to_deseq2(ps_fungi_clean_chao1_12, ~ Chao1_dif_cat)

# Convert counts to integer
vst_phyloseq_dseq_fungi = estimateSizeFactors(vst_phyloseq_dseq_fungi, geoMeans = apply(counts(vst_phyloseq_dseq_fungi), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(vst_phyloseq_dseq_fungi, blind = TRUE)
vst_blind_mat_fungi <- SummarizedExperiment::assay(vst_blind) 
vst_blind_mat_fungi <- t(vst_blind_mat_fungi) 
vst_blind_mat_fungi[which(vst_blind_mat_fungi < 0)] <- 0 
dists <- dist(t(assay(vst_phyloseq_dseq_fungi)))

# Computing Bray-Curtis Dissimilarities and PCoA
comm_vst_blind_mat_fungi<- vegdist(vst_blind_mat_fungi, "bray")
PCoA_comm_vst_blind_mat_fungi <- capscale(comm_vst_blind_mat_fungi ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat_fungi$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat_fungi$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat_fungi)$sites

# Filter dataframe to remove unchanged Chao1 values
fungi_data_12_chao1_clean <- fungi_data_12_chao1 %>% filter(Chao1_dif_cat != "Unchanged")

# Save scores into metadata tables
row.names(fungi_data_12_chao1_clean) == row.names(scores(PCoA_comm_vst_blind_mat_fungi)$sites)
fungi_data_12_chao1_clean$PCoA1_fungi_12 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,1]
fungi_data_12_chao1_clean$PCoA2_fungi_12 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,2]

# Relabel column values for Chao1 trend beta-diversity plot
fungi_data_12_chao1_clean$Chao1_dif_cat_num <- fungi_data_12_chao1_clean$Chao1_dif_cat
fungi_data_12_chao1_clean <- fungi_data_12_chao1_clean %>% mutate(Chao1_dif_cat_num = dplyr::recode(Chao1_dif_cat_num,
                                                                                      "Increase" = "Increase (n = 25)",
                                                                                      "Decrease" = "Decrease (n = 63)"))

# PCoA plot for beta diversity at 12 months by Chao1 trend (Figure 1D)
PCoA_fungi_chao1_12 <- qplot(PCoA1_fungi_12, PCoA2_fungi_12, xlab = "PCoA1", ylab = "PCoA2",
                    size = 0.5, shape = Sample_time, fill = Chao1_dif_cat, color = Chao1_dif_cat, data = (fungi_data_12_chao1_clean))

PCoA_fungi_chao1_12 + stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Chao1_dif_cat, color = Chao1_dif_cat)) +
  scale_shape_manual(name = "Age (Months)", values = c(22:24)) + 
  theme_bw() + 
  xlab("PCoA1 (13.9%)")+
  ylab("PCoA2 (8.5%)")+
  ggtitle("12 Months")+
  ylim(-2,2)+
  xlim(-2,3)+
  theme(legend.title = element_text(colour = "black", size = 9, face = "bold"),
        legend.text = element_text(colour = "black", size = 9),
        legend.position = "none",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Richness Pattern", values = c("black","black")) + 
  scale_fill_manual(name = "Richness Pattern", values = c("#2db27d", "#481b6d")) + 
  guides(size = FALSE) 

# PERMANOVA by Chao1 trend at 12 months
permanova_fungi <- adonis(comm_vst_blind_mat_fungi ~ Chao1_dif_cat, data = fungi_data_12_chao1_clean, permutations = 999)
permanova_fungi$aov.tab

# Repeat beta-diversity for Shannon trend (Figure S1B)

# Filter ps object by sample time and to remove NA values for Shannon
ps_fungi_clean_shannon <- subset_samples(ps_fungi_clean, Shannon_dif_cat != "NA")
ps_fungi_clean_shannon_3 <- subset_samples(ps_fungi_clean_shannon, Sample_time == 3)

# Variance Stabilizing Transformation
vst_phyloseq_dseq_fungi <- phyloseq_to_deseq2(ps_fungi_clean_shannon_3, ~Shannon_dif_cat)

# Convert counts to integer
vst_phyloseq_dseq_fungi = estimateSizeFactors(vst_phyloseq_dseq_fungi, geoMeans = apply(counts(vst_phyloseq_dseq_fungi), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(vst_phyloseq_dseq_fungi, blind = TRUE)
vst_blind_mat_fungi <- SummarizedExperiment::assay(vst_blind) # Extract transformed OTU table
vst_blind_mat_fungi <- t(vst_blind_mat_fungi) # Transpose data
vst_blind_mat_fungi[which(vst_blind_mat_fungi < 0)] <- 0 # Set counts less than 0 to zero
dists <- dist(t(assay(vst_phyloseq_dseq_fungi)))

# Computing Bray-Curtis Dissimilarities and PCoA
comm_vst_blind_mat_fungi<- vegdist(vst_blind_mat_fungi, "bray")
PCoA_comm_vst_blind_mat_fungi <- capscale(comm_vst_blind_mat_fungi ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat_fungi$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat_fungi$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat_fungi)$sites

# Save scores into metadata tables
row.names(fungi_data_3_shannon) == row.names(scores(PCoA_comm_vst_blind_mat_fungi)$sites)
fungi_data_3_shannon$PCoA1_fungi_3 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,1]
fungi_data_3_shannon$PCoA2_fungi_3 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,2]

# Relabel column values for Shannon trend beta-diversity plot
fungi_data_3_shannon$Shannon_dif_cat_num <- fungi_data_3_shannon$Shannon_dif_cat
fungi_data_3_shannon <- fungi_data_3_shannon %>% mutate(Shannon_dif_cat_num = dplyr::recode(Shannon_dif_cat_num,
                                                                                            "Increase" = "Increase (n = 18)",
                                                                                            "Decrease" = "Decrease (n = 73)"))

# PCoA plot for beta diversity at 3 months by Shannon trend (Figure S1B)
PCoA_fungi_shannon_3 <- qplot(PCoA1_fungi_3, PCoA2_fungi_3, xlab = "PCoA1", ylab = "PCoA2",
                    size = 0.5, shape = Sample_time, fill = Shannon_dif_cat, color = Shannon_dif_cat, data = (fungi_data_3_shannon))

PCoA_fungi_shannon_3 + stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Shannon_dif_cat, color = Shannon_dif_cat)) +
  scale_shape_manual(name = "Age (Months)", values = 21) + #shape selection
  theme_bw() + 
  xlab("PCoA1 (6.2%)")+
  ylab("PCoA2 (5.0%)")+
  ggtitle("3 Months")+
  #ylim(-2,2)+
  #xlim(-2,3)+
  theme(legend.title = element_text(colour = "black", size = 9, face = "bold"),
        legend.text = element_text(colour = "black", size = 9),
        legend.position = "bottom",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Alpha Diveristy Pattern", values = c("black","black")) + #outer shape color
  scale_fill_manual(name = "Alpha Diveristy Pattern", values = c("#2db27d", "#481b6d")) + #inside shape color
  guides(size = FALSE) # To hide legend, add this inside parenthesis: color = FALSE, fill = FALSE, shape = FALSE

# PERMANOVA by Shannon trend at 3 months
permanova_fungi <- adonis(comm_vst_blind_mat_fungi ~ Shannon_dif_cat, data = fungi_data_3_shannon, permutations = 999)
permanova_fungi$aov.tab

# Repeat for Shannon trend at 12 months

# Filter ps object by sample time and to remove NA values for Shannon
ps_fungi_clean_shannon_12 <- subset_samples(ps_fungi_clean_shannon, Sample_time == 12)

# Variance Stabilizing Transformation
vst_phyloseq_dseq_fungi <- phyloseq_to_deseq2(ps_fungi_clean_shannon_12, ~Shannon_dif_cat)

# Convert counts to integer
vst_phyloseq_dseq_fungi = estimateSizeFactors(vst_phyloseq_dseq_fungi, geoMeans = apply(counts(vst_phyloseq_dseq_fungi), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(vst_phyloseq_dseq_fungi, blind = TRUE)
vst_blind_mat_fungi <- SummarizedExperiment::assay(vst_blind) # Extract transformed OTU table
vst_blind_mat_fungi <- t(vst_blind_mat_fungi) # Transpose data
vst_blind_mat_fungi[which(vst_blind_mat_fungi < 0)] <- 0 # Set counts less than 0 to zero
dists <- dist(t(assay(vst_phyloseq_dseq_fungi)))

# Computing Bray-Curtis Dissimilarities and PCoA
comm_vst_blind_mat_fungi<- vegdist(vst_blind_mat_fungi, "bray")
PCoA_comm_vst_blind_mat_fungi <- capscale(comm_vst_blind_mat_fungi ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat_fungi$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat_fungi$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat_fungi)$sites

# Save scores into metadata table
row.names(fungi_data_12_shannon) == row.names(scores(PCoA_comm_vst_blind_mat_fungi)$sites)
fungi_data_12_shannon$PCoA1_fungi_12 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,1]
fungi_data_12_shannon$PCoA2_fungi_12 <- scores(PCoA_comm_vst_blind_mat_fungi)$sites[,2]

# Relabel column values for Shannon trend beta-diversity plot
fungi_data_12_shannon$Shannon_dif_cat_num <- fungi_data_12_shannon$Shannon_dif_cat
fungi_data_12_shannon <- fungi_data_12_shannon %>% mutate(Shannon_dif_cat_num = dplyr::recode(Shannon_dif_cat_num,
                                                                                              "Increase" = "Increase (n = 18)",
                                                                                              "Decrease" = "Decrease (n = 73)"))

# PCoA plot for beta diversity at 12 months by Shannon trend (Figure S1B)
PCoA_fungi_shannon_12 <- qplot(PCoA1_fungi_12, PCoA2_fungi_12, xlab = "PCoA1", ylab = "PCoA2",
                    size = 0.5, shape = Sample_time, fill = Shannon_dif_cat_num, color = Shannon_dif_cat_num, data = (fungi_data_12_shannon))

PCoA_fungi_shannon_12 + stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Shannon_dif_cat_num, color = Shannon_dif_cat_num)) +
  scale_shape_manual(name = "Age (Months)", values = c(22)) + 
  theme_bw() + 
  xlab("PCoA1 (13.5%)")+
  ylab("PCoA2 (8.6%)")+
  ggtitle("12 Months")+
  ylim(-2,3)+
  xlim(-1.75,2.75)+
  theme(legend.title = element_text(colour = "black", size = 9, face = "bold"),
        legend.text = element_text(colour = "black", size = 9),
        legend.position = "none",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 10, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Alpha Diveristy Pattern", values = c("black","black")) +
  scale_fill_manual(name = "Alpha Diveristy Pattern", values = c("#2db27d", "#481b6d")) + 
  guides(size = FALSE) 

# PERMANOVA by Shannon trend at 12 months
permanova_fungi <- adonis(comm_vst_blind_mat_fungi ~ Shannon_dif_cat, data = fungi_data_12_shannon, permutations = 999)
permanova_fungi$aov.tab

## Relative Abundance of Top 10 Fungal Genera (Figures 1E-F & S1C)

# Create ps object and data frame of top ten genera for plotting relative abundance
ps_fungi_genus <- tax_glom(ps_fungi_clean, taxrank = "Genus")
ps_fungi_genus <- transform_sample_counts(ps_fungi_genus, function(OTU) OTU/sum(OTU))
pmelt_fungi_genus <- psmelt(ps_fungi_genus)
sum(pmelt_fungi_genus$Abundance)

top10genus_fungi <- names(sort(taxa_sums(ps_fungi_genus), decreasing = TRUE))[1:10]
ps_top10genus_fungi <- prune_taxa(top10genus_fungi, ps_fungi_genus)
ps_top10genus_fungi
pmelt_fungi_genus10 <- psmelt(ps_top10genus_fungi)
unique(pmelt_fungi_genus10$Genus)

# Determine percent of community top 10 genera account for at 3 and 12 months (mean and standard deviation)
ps_fungi_genus_3_top10 <- subset_samples(ps_top10genus_fungi, Sample_time == 3)
asv_genus_3 = as.data.frame(otu_table(ps_fungi_genus_3_top10), "matrix") %>% mutate(Other = 100*(1-rowSums(.)))
100-mean(asv_genus_3$Other)
sd(asv_genus_3$Other)

ps_fungi_genus_12_top10 <- subset_samples(ps_top10genus_fungi, Sample_time == 12)
asv_genus_12 = as.data.frame(otu_table(ps_fungi_genus_12_top10), "matrix") %>% mutate(Other = 100*(1-rowSums(.)))
100-mean(asv_genus_12$Other)
sd(asv_genus_12$Other)

# Rename genera for tidy plots by removing "g__" values from names
pmelt_fungi_genus10$Genus_new <- pmelt_fungi_genus10$Genus
pmelt_fungi_genus10 <- pmelt_fungi_genus10 %>% mutate(Genus_new = str_sub(Genus_new, 4, -1))
pmelt_fungi_genus10$Genus <- pmelt_fungi_genus10$Genus_new 

# Plot relative abundance of top 10 fungal genera by sample time (Figure 1E)
fungal_genus <- ggplot(pmelt_fungi_genus10, aes(x = Sample_time, y = Abundance, fill = Genus, colour = Genus)) + 
  geom_bar(stat = "identity", position = "fill") + 
  xlab("Age (Months)") + 
  ylab("Relative Abundance") +
  scale_fill_manual(values=c("#440154","mediumorchid4","#0093af","deepskyblue4","darkcyan","#35b779","olivedrab3","#cde11d","gold1","#fde725"))+
  scale_color_manual(values = c("#440154","mediumorchid4","#0093af","deepskyblue4","darkcyan","#35b779","olivedrab3","#cde11d","gold1","#fde725")) +
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9),
        legend.text.align = 0,
        legend.position = "none",
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))
fungal_genus

# Repeat for relative abundance of top 10 genera by Chao1 trend

# Relabel sample time columns
pmelt_fungi_genus10$Sample_time_month <- pmelt_fungi_genus10$Sample_time
pmelt_fungi_genus10 <- pmelt_fungi_genus10 %>% mutate(Sample_time_month = dplyr::recode(Sample_time_month,
                                                                                        "3" = "3 Months",
                                                                                        "12" = "12 Months"))

# Filter out NA and unchanged rows for Chao1 trend & reorder categories
pmelt_fungi_genus10_chao1 <- pmelt_fungi_genus10 %>% filter(Chao1_dif_cat == c("Increase", "Decrease"))
pmelt_fungi_genus10_chao1$Chao1_dif_cat <- factor(pmelt_fungi_genus10_chao1$Chao1_dif_cat, levels = c("Decrease", "Increase"))

# Plot relative abundance of top 10 fungal genera by sample time and Chao1 trend (Figure 1F)
fungal_genus_chao1 <- ggplot(pmelt_fungi_genus10_chao1, aes(x = Chao1_dif_cat, y = Abundance, fill = Genus, colour = Genus)) + 
  geom_bar(stat = "identity", position = "fill") + 
  xlab("Richness Pattern") + 
  ylab("Relative Abundance") +
  facet_wrap(~ Sample_time_month) +
  scale_fill_manual(values=c("#440154","mediumorchid4","#0093af","deepskyblue4","darkcyan","#35b779","olivedrab3","#cde11d","gold1","#fde725"))+
  scale_color_manual(values = c("#440154","mediumorchid4","#0093af","deepskyblue4","darkcyan","#35b779","olivedrab3","#cde11d","gold1","#fde725")) +
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9, face = "italic"),
        legend.text.align = 0,
        legend.position = "right",
        strip.text.x = element_text(size = 9,face = "bold"),
        strip.background = element_rect(color = "white", fill = "white", size = 0.5),
        plot.title = element_text(hjust = 0.5))
fungal_genus_chao1

# Repeat for relative abundance of top 10 genera by Shannon trend

# Filter out NA rows for Shannon trend & reorder categories
pmelt_fungi_genus10_shannon <- pmelt_fungi_genus10 %>% filter(Shannon_dif_cat == c("Increase", "Decrease"))
pmelt_fungi_genus10_shannon$Shannon_dif_cat <- factor(pmelt_fungi_genus10_shannon$Shannon_dif_cat, levels = c("Decrease", "Increase"))

# Plot relative abundance of top 10 fungal genera by sample time and Shannon trend (Figure S1C)
fungal_genus_shannon <- ggplot(pmelt_fungi_genus10_clean, aes(x = Shannon_dif_cat, y = Abundance, fill = Genus, colour = Genus)) + 
  geom_bar(stat = "identity", position = "fill") + 
  xlab("Alpha Diversity Pattern") + 
  ylab("Relative Abundance") +
  facet_wrap(~ Sample_time_month) +
  scale_fill_manual(values=c("#440154","mediumorchid4","#0093af","deepskyblue4","darkcyan","#35b779","olivedrab3","#cde11d","gold1","#fde725"))+
  scale_color_manual(values = c("#440154","mediumorchid4","#0093af","deepskyblue4","darkcyan","#35b779","olivedrab3","#cde11d","gold1","#fde725")) +
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 9),
        legend.title = element_text(face = "bold", size = 9), 
        legend.text = element_text(size = 9, face = "italic"),
        legend.text.align = 0,
        legend.position = "right",
        strip.text.x = element_text(size = 9,face = "bold"),
        strip.background = element_rect(color = "white", fill = "white", size = 0.5),
        plot.title = element_text(hjust = 0.5))
fungal_genus_shannon

# Relative abundance plots per genera by sample time and Chao1 pattern for top 10 genera (Figure S2)

# Filter genus-level data per genera for individual boxplots of relative abundance by Chao1 pattern at 3 and 12 months
pmelt_candida <- pmelt_fungi_genus_chao1 %>% filter(Genus == "g__Candida")

# Plot relative abundance per genera by Chao1 pattern at 3 and 12 months (Figure S2)
candida_fig <- ggplot(pmelt_candida, aes(x = Chao1_dif_cat, y = Abundance, fill = Chao1_dif_cat, color = Chao1_dif_cat)) + 
  geom_boxplot(alpha = 0.6) +
  geom_jitter(position = position_jitter(0.2),  size = 1.2) + 
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("#2db27d", "#481b6d"))+
  xlab("Richness Pattern")+
  ylab("Relative Abundance")+
  ggtitle("Candida spp.")+
  facet_wrap(~Sample_time_month)+
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  theme(axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(size = 9),
        strip.text.x = element_text(size = 10,face = "bold"),
        strip.background = element_rect(color = "white", fill = "white", size = 0.5),
        plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 10))
candida_fig

# Load clr-transformed abundance ps object previously generated for taxa heat maps
load("ps_fungi_clean_clr")

# Remove participants with NA or Unchanged values for Chao1 pattern from ps object
ps_fungi_genus_noNA <- subset_samples(ps_fungi_clean_clr, Chao1_dif_cat != "Unchanged")
ps_fungi_genus_noNA # Check N
ps_fungi_genus_noNA <- subset_samples(ps_fungi_genus_noNA, Chao1_dif_cat != "NA")
ps_fungi_genus_noNA # Check N

# Create dataframe with abundance information per genera for top 10 fungal genera by sample time
ps_candida <- subset_taxa(ps_fungi_genus_noNA, Genus == "g__Candida")
pmelt_candida <- psmelt(ps_candida)
pmelt_candida_3 <- pmelt_candida %>% filter(Sample_time == 3)
pmelt_candida_12 <- pmelt_candida %>% filter(Sample_time == 12)

# Perform Shapiro test for normality and Wilcox test per sample time for differences in abundance by Chao1 pattern
shapiro.test(pmelt_candida_3$Abundance)
wilcox.test(pmelt_candida_3$Abundance ~ Chao1_dif_cat, data = pmelt_candida_3)
shapiro.test(pmelt_candida_12$Abundance)
wilcox.test(pmelt_candida_12$Abundance ~ Chao1_dif_cat, data = pmelt_candida_12)

# Repeat for Cladosporium, Debaromyces, Fomitopsis, Ganoderma, Malassezia, Mycosphaerella, Resinicium, Rhodotorula & Saccharomyces
