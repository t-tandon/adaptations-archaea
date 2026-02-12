#import file
library(readr)
df <- read_tsv("../feature_set.tsv")
#numeric only
nm_df <- df[,-1]

#==SCATTERPLOT OF FEATURES== 
#load library
library(tidyr)
library(ggplot2)

#reshape into long format
df_long <- pivot_longer(nm_df, cols = -ogt, names_to = "trait", values_to = "value")

#ggplot with facetwrap
ggplot(df_long, aes(x = ogt, y = value)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~trait, scales = "free_y") +
  theme_minimal()

#==CORRELATION MATRIX== 
cor_matrix <- cor(nm_df, use="pairwise.complete.obs")
round(cor_matrix, 2)
heatmap(cor_matrix, symm=TRUE)

#RUN SIMPLE LINEAR REGRESSION 
model_lm <- lm(ogt ~ ., data = nm_df)
summary(model_lm)

#RUN PHYLOGENETIC CORRELATION ANALYSIS
library(ape)
library(nlme)
library(caper)
library(phytools)

#load tree
tree <- read.tree("archaeal_tree")
#remove all node labels
tree$node.label <- NULL
#midpoint rooting
tree_root <- midpoint.root(tree)
#load data
data <- read.csv("feature_set.tsv",sep="\t",row.names = 1)
#clean tree labels
tree_root$tip.label <- sub(".*((GCA|GCF)_\\d+(\\.\\d+)?).*", "\\1", tree_root$tip.label)
tree_root$tip.label <- sub("\\.\\d+$", "", tree_root$tip.label)

#sync tree and table
common <- intersect(tree_root$tip.label, rownames(data))
tree_root <- drop.tip(tree_root, setdiff(tree_root$tip.label, common))
data <- data[common, ]

data$species <- rownames(data) 
comp_data <- comparative.data(tree_root, data, names.col = "species", vcv = TRUE, na.omit = FALSE)
#formula
formula <- ogt ~ frac_ilvwyg + ratio_charged_polar + trna_gc + aromo + nc + genome_size + gene_density + gc3s

#run PGLS
pgls_model <- pgls(formula, data = comp_data, lambda = "ML") 
summary(pgls_model)
#Biological conclusion robust: amino acid composition + trna gc content are consistent predictors
