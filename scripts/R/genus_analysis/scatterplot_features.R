#load libraries
library(ape)
library(nlme)
library(dplyr)
library(ggplot2)
library(caper)
library(phytools)

#load data
data <- read.csv("../feature_set_phylum_class.tsv",sep="\t",row.names = 1)

features <- c("frac_ilvwygerkp", "ratio_charged_polar", "trna_gc", 
  "aromo", "nc", "genome_size", "gene_density", 
  "gc3s", "frac_hydrophobic")

#convert phylum + class to factors
data$phylum <- as.factor(data$phylum)
data$class  <- as.factor(data$class)

#==PHYLUM== 
#only for phylum > 20 likely to be statistically significant
phyla_counts <- table(data$phylum)
good_phyla <- names(phyla_counts[phyla_counts >= 20])

#subset data
data_filtered <- subset(data, phylum %in% good_phyla)

#long format
long_dat <- data_filtered %>%
  pivot_longer(
    cols = all_of(features),
    names_to = "feature",
    values_to = "value")

#exploratory figure
p <- ggplot(long_dat, aes(x = value, y = ogt, color = phylum)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~ feature, scales = "free_x", ncol = 3) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(size = 12, face = "bold"), legend.position = "bottom") +
  labs(
    x = "Feature value",
    y = "Optimal growth temperature (°C)",
    color = "Phylum",
    title = "Exploratory plots: OGT vs genomic features"
  )

ggsave("ogt_vs_features_phylum_filtered.pdf", p, width = 20, height = 30)

##==CLASS== 
class_counts <- table(data$class)
good_class <- names(class_counts[class_counts >= 10 ])

#subset data
data_filtered <- subset(data, class %in% good_class)

#long format
long_dat <- data_filtered %>%
  pivot_longer(
    cols = all_of(features),
    names_to = "feature",
    values_to = "value")

#exploratory figure
p <- ggplot(long_dat, aes(x = value, y = ogt, color = class)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  facet_wrap(~ feature, scales = "free_x", ncol = 3) +
  theme_bw(base_size = 13) +
  theme(strip.text = element_text(size = 12, face = "bold"), legend.position = "bottom") +
  labs(
    x = "Feature value",
    y = "Optimal growth temperature (°C)",
    color = "Class",
    title = "Exploratory plots: OGT vs genomic features"
  )

ggsave("class_ogt_vs_features-filtered.pdf", p, width = 20, height = 30)
