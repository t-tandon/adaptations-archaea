#LOAD FILES 
merged <- read.csv("merged_summary.tsv", sep="\t")

## RUN GLS ANALYSIS
library(ape)
library(nlme)
library(dplyr)

# Load tree
tree <- read.tree("archaeal_tree")
tree$node.label <- NULL
tree <- midpoint.root(tree)

# Clean accession labels
tree$tip.label <- sub(".*((GCA|GCF)_\\d+(\\.\\d+)?).*", "\\1", tree$tip.label)
tree$tip.label <- sub("\\.[0-9]+$", "", tree$tip.label)

#create dataframe
data_min <- merged[, c("accession", "class", "ogt", "genome_size",
                       "gc_content", "frac_ilvwygerkp",
                       "ratio_charged_polar", "trna_gc")]
na.omit(data_min)


#match tree/data
common <- intersect(tree$tip.label, data_min$accession)
tree2 <- drop.tip(tree, setdiff(tree$tip.label, common))
data2 <- data_min %>% filter(accession %in% common)
# reorder rows to match tree tip order
data2 <- data2[match(tree2$tip.label, data2$accession), ]

# sanity check
all(tree2$tip.label == data2$accession)

#identify classes with enough power
min_n <- 20
class_counts <- table(data2$class)
good_class <- names(class_counts[class_counts >= min_n])
good_class

## RUN GLS PER CLASS
results_table <- data.frame()

for (cl in good_class) {
  
  cat("Running:", phy, "\n")
  
  sub_data <- data2 %>% filter(class == cl)
  
  sub_tree <- drop.tip(tree2,
                       setdiff(tree2$tip.label,
                               sub_data$accession))
  
  sub_data <- sub_data[match(sub_tree$tip.label,
                             sub_data$accession), ]
  
  if (nrow(sub_data) < min_n) next
  
  model <- try(
    gls(
      ogt ~ genome_size + gc_content +
        frac_ilvwygerkp +
        ratio_charged_polar +
        trna_gc,
      data = sub_data,
      correlation = corPagel(1, sub_tree, fixed = FALSE),
      method = "ML"
    ),
    silent = TRUE
  )
  
  if (inherits(model, "try-error")) next
  
  summ <- summary(model)
  
  lambda_est <- coef(model$modelStruct$corStruct, unconstrained = FALSE)
  
  coef_table <- as.data.frame(summ$tTable)
  coef_table$term <- rownames(coef_table)
  coef_table$phylum <- phy
  coef_table$n <- nrow(sub_data)
  coef_table$lambda <- lambda_est
  coef_table$AIC <- AIC(model)
  
  results_table <- rbind(results_table, coef_table)
}

# SAVE SUMMARY TABLE 
write.table(results_table, "per_class_analysis.tsv", sep="\t",
            row.names=FALSE, quote=FALSE)

