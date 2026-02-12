#load libraries
library(ape) 
library(geiger) 
library(tidyverse)
library(phytools) 



#load tree and data
tree <- read.tree("archaeal_tree")
#clean tip labels
tree$tip.label <- sub(".*((GCA|GCF)_\\d+(\\.\\d+)?).*", "\\1", tree$tip.label)
tree$tip.label <- sub("\\.\\d+$", "", tree$tip.label)
#root tree
tree <- midpoint.root(tree)
tree$node.label <- NULL
#load metadata, match to tree
data <- read.table("feature_set_phylum_class.tsv",
                   sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
#keep only common 
common <- intersect(tree$tip.label, rownames(data))
tree2 <- drop.tip(tree, setdiff(tree$tip.label, common))
data2 <- data[tree$tip.label, , drop = FALSE]
#replace with new labels
old_labels <- tree2$tip.label
new_labels <- paste0(data2[old_labels,"phylum"], " (",
                     data2[old_labels, "class"], ") -", 
                     old_labels)
#apply to tree
tree2$tip.label <- new_labels
rownames(data2) <- new_labels

#traits to analyse
test_traits <- c("ogt", "genome_size", "gene_density", "nc", "ratio_charged_polar","frac_ilvwygerkp", 
                 "frac_hydrophobic", "trna_gc", "gc3s", "aromo")

#prepare named vector
tip_values <- list()

for (trait in test_traits) {
  v <- data2[[trait]]
  names(v) <- rownames(data2)
  v <- v[!is.na(v)]     # geiger breaks with NA values
  tip_values[[trait]] <- v
}

##MODEL TESTING
continuous_models <- c("BM","EB","lambda","kappa","delta","white")

geiger_results <- tibble()

#fit all traits across all models

for (trait in names(tip_values)) {
  cat("\n===== Trait:", trait, "=====\n")
  used_trait <- tip_values[[trait]]
  # tree must match taxon set of trait
  t <- drop.tip(tree2, setdiff(tree2$tip.label, names(used_trait)))
  for (model in continuous_models) {
    cat("Model:", model, "\n")
    fit <- fitContinuous(
      phy = t,
      dat = used_trait,
      model = model,
      control = list(
        method = c("subplex", "L-BFGS-B", "Nelder-Mead"),
        niter = 500,
        FAIL = 1e+200,
        hessian = FALSE
      )
    )
    stats <- fit$opt %>%
      as_tibble() %>%
      mutate(
        trait = trait,
        model = model
      )
    
    geiger_results <- bind_rows(geiger_results, stats)
  }} 

best_models <- geiger_results %>%
  group_by(trait) %>%
  slice_min(aicc, n = 1) %>%
  ungroup()

# Parameter extractor
get_param <- function(row) {
  model <- row$model
  
  if (model == "lambda") return(row$lambda)
  if (model == "kappa")  return(row$kappa)
  if (model == "delta")  return(row$delta)
  if (model == "EB")     return(row$a)
  return(NA_real_)        # BM, white, no parameter
}

best_models <- best_models %>%
  rowwise() %>%
  mutate(param = get_param(cur_data())) %>%
  ungroup()

#save reuslts
write_tsv(best_models, "best_model_per_trait_archaea.txt")

#create output folders
dir.create("ASR_plots", showWarnings = FALSE)

#storage objects
ancestral_states <- list() 
rescaled_trees <- list() 
contmaps <- list() 

#loop over all traits
for (trait_name in names(tip_values)) {
  
  message("==Processing Trait for", trait_name, "==")
  
  trait_vector <- tip_values[[trait_name]]
  best_model <- best_mods %>% filter(trait == trait_name)
  
  #rescale to model
  if (best_model$model %in% c("lambda", "kappa", "delta")){
    parameter_value <- best_model[[best_model$model]]
    t_rescaled <- phytools::rescale(tree2, 
                                    model = best_model$model, parameter_value)
  } else if (best_model$model == "EB") {
    t_rescaled <- phytools::rescale(tree2, 
                                    model = "EB", 
                                    best_model$a)
  } else if (best_model$model == "BM") {
    t_rescaled <- tree2
  }
  rescaled_trees[[trait_name]] <- t_rescaled
  
  #ancestral reconstruction
  anc <- anc.recon(tree = t_rescaled,
                   trait_data = trait_vector, vars = TRUE, CI = TRUE)
  
  ancestral_states[[trait_name]] <- anc
  
  #visualisation with contmap
  cmap <- phytools::contMap(t_rescaled, trait_vector, plot = FALSE)
  #save plot
  pdf(file = paste0("ASR_plots/", trait_name, "_contMap.pdf"), 
      width = 10, height = 50)
  plot(cmap, 
       fsize = 0.3, 
       lwd = 4, 
       outline = FALSE, 
       legend = TRUE)
  dev.off() 
  contmaps[[trait_name]] <- cmap
}

saveRDS(ancestral_states, "ASR_archaea.rds")
saveRDS(rescaled_trees, "ASR_trees_rescaled.rds")
saveRDS(contmaps, "ASR_contmaps.rds")
message("==Pipeline complete!")



