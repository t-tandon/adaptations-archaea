#load libraries
library(ape)
library(nlme)
library(dplyr)
library(ggplot2)
library(caper)
library(phytools)

#PROCESS TREE

#load tree
tree <- read.tree("archaeal_tree")
# remove all node labels
tree$node.label <- NULL
#midpoint rooting
tree_root <- midpoint.root(tree)
#load data
data <- read.csv("feature_set_phylum_class.tsv",sep="\t",row.names = 1)
#clean tree labels
tree_root$tip.label <- sub(".*((GCA|GCF)_\\d+(\\.\\d+)?).*", "\\1", tree_root$tip.label)
tree_root$tip.label <- sub("\\.\\d+$", "", tree_root$tip.label)
#add species column 
data$species <- rownames(data)
#define formula
formula <- ogt ~ frac_ilvwygerkp + ratio_charged_polar + trna_gc + aromo + nc + genome_size + gene_density + gc3s + frac_hydrophobic

#PER PHYLUM PGLS
#list of phyla with > 20 sequences 
phylum_counts <- table(data$phylum)
all_phyla <- names(phylum_counts[phylum_counts >= 20])

#to store results
results_list = list() 

#loop thru phyla
for (phy in all_phyla) {
  cat("Running PGLS for ", phy, "...")
  #subset data
  sub_data <- subset(data, phylum == phy)
  #prune
  sub_tree <- drop.tip(tree_root, setdiff(tree_root$tip.label, sub_data$species))
  #built comparative data obj 
  sub_comp <- comparative.data(sub_tree, sub_data, 
                               names.col = "species", 
                               vcv = TRUE, na.omit = FALSE)
  #fit pgls
  sub_model <- pgls(formula, data = sub_comp, lambda = "ML")
  #store summary
  results_list[[phy]] <- summary(sub_model)
}

print(results_list)

##CLASS PGLS
#get all classes with at least N >= 10
class_counts <- table(data$class)
all_classes <- names(class_counts[class_counts >= 10])

#PGLS per class
class_results <- list() 
for (cl in all_classes) {
  cat("Running PGLS for ", cl, "...")
  
  #subset data
  sub_data <- subset(data, class == cl)
  #prune
  sub_tree <- drop.tip(tree_root, setdiff(tree_root$tip.label, sub_data$species))
  #built comparative data obj 
  sub_comp <- comparative.data(sub_tree, sub_data, 
                               names.col = "species", 
                               vcv = TRUE, na.omit = FALSE)
  #fit pgls
  sub_model2 <- pgls(formula, data = sub_comp, lambda = "ML")
  #store summary
  class_results[[cl]] <- summary(sub_model2)
}
print(class_results)
