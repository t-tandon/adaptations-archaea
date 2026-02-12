#import file
library(readr)
df <- read_tsv("summary.tsv")

#numeric only
nm_df <- df[,-1]

## INITIAL EXPLORATORY ANALYSIS

#scatterplot
library(tidyr)

#reshape into long format
df_long <- pivot_longer(nm_df, cols = -ogt, names_to = "trait", values_to = "value")

library(ggplot2)

#ggplot with facetwrap
ggplot(df_long, aes(x = ogt, y = value)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~trait, scales = "free_y") +
  theme_minimal()

#correlation matrix
cor_matrix <- cor(nm_df, use="pairwise.complete.obs")
round(cor_matrix, 2)
#heatmap
heatmap(cor_matrix, symm=TRUE)


#correlation coefficients
library(car)
vif(lm(ogt ~ ., data = nm_df))

## SIMPLE GLS MODEL
#load libraries
library(ape)
library(caper)
library(phytools)

#load tree
tree <- read.tree("archaeal_tree")
tree$node.label <- NULL
tree <- midpoint.root(tree)

#clean labels
tree$tip.label <- sub(".*((GCA|GCF)_\\d+(\\.\\d+)?).*", "\\1", tree$tip.label)
tree$tip.label <- sub("\\.\\d+$", "", tree$tip.label)

#load data
data <- read.csv("summary.tsv", sep="\t", row.names=1)

#keep only needed predictors
data_min <- data[, c("ogt", "genome_size",
                     "gc_content", "frac_ilvwygerkp",
                     "ratio_charged_polar", "trna_gc")]

data_min <- na.omit(data_min)

#match tree/data
common <- intersect(tree$tip.label, rownames(data_min))
tree <- drop.tip(tree, setdiff(tree$tip.label, common))
data_min <- data_min[common, ]

data_min$species <- rownames(data_min)

library(nlme)

gls_model <- gls(
  ogt ~ genome_size + gc_content + frac_ilvwygerkp + ratio_charged_polar + trna_gc,
  data = data_min,
  correlation = corPagel(1, tree, fixed=FALSE),
  method="ML"
)

summary(gls_model)

write.table(data_min,
            "summary_clean",
            sep="\t",
            quote=FALSE,
            row.names=FALSE)
