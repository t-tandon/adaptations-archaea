# =========================
# LOAD LIBRARIES
# =========================

library(ape)
library(phytools)
library(geiger)
library(dplyr)
library(tibble)
library(readr)

# =========================
# LOAD AND CLEAN TREE
# =========================

tree <- read.tree("archaeal_tree")

# clean accession labels
tree$tip.label <- sub(".*((GCA|GCF)_\\d+(\\.\\d+)?).*", "\\1", tree$tip.label)
tree$tip.label <- sub("\\.[0-9]+$", "", tree$tip.label)

tree <- midpoint.root(tree)
tree$node.label <- NULL

# =========================
# LOAD MERGED SPECIES DATA
# =========================

data <- read.delim("merged_summary.tsv",
                   sep = "\t",
                   header = TRUE,
                   stringsAsFactors = FALSE)

# =========================
# DEFINE TRAITS
# =========================

test_traits <- c("ogt",
                 "genome_size",
                 "ratio_charged_polar",
                 "frac_ilvwygerkp", "gc_content", 
                 "trna_gc")

# =========================
# CLEAN DATA
# =========================

data_clean <- data %>%
  select(accession, phylum, class, all_of(test_traits)) %>%
  filter(!is.na(accession))

# =========================
# MATCH TREE AND DATA
# =========================

common <- intersect(tree$tip.label, data_clean$accession)

tree2 <- drop.tip(tree, setdiff(tree$tip.label, common))

data2 <- data_clean %>%
  filter(accession %in% common)

# reorder to match tree order
data2 <- data2[match(tree2$tip.label, data2$accession), ]

stopifnot(all(tree2$tip.label == data2$accession))

# =========================
# RENAME TREE TIPS WITH TAXONOMY
# =========================

new_labels <- paste0(data2$phylum,
                     " (", data2$class, ") - ",
                     data2$accession)

tree2$tip.label <- new_labels
rownames(data2) <- new_labels

# =========================
# PREPARE TRAIT VECTORS
# =========================

tip_values <- list()

for (trait in test_traits) {
  v <- data2[[trait]]
  names(v) <- rownames(data2)
  v <- v[!is.na(v)]
  tip_values[[trait]] <- v
}

# =========================
# FIT EVOLUTIONARY MODELS
# =========================

continuous_models <- c("BM","EB","lambda","kappa","delta","white")

geiger_results <- tibble()

for (trait in names(tip_values)) {
  
  cat("\n===== Trait:", trait, "=====\n")
  
  used_trait <- tip_values[[trait]]
  
  t_sub <- drop.tip(tree2,
                    setdiff(tree2$tip.label,
                            names(used_trait)))
  
  for (model in continuous_models) {
    
    cat("Model:", model, "\n")
    
    fit <- fitContinuous(
      phy = t_sub,
      dat = used_trait,
      model = model,
      control = list(
        method = c("subplex", "L-BFGS-B", "Nelder-Mead"),
        niter = 500,
        FAIL = 1e+200,
        hessian = FALSE
      )
    )
    
    stats <- as_tibble(fit$opt) %>%
      mutate(trait = trait,
             model = model)
    
    geiger_results <- bind_rows(geiger_results, stats)
  }
}

# =========================
# SAVE MODEL RESULTS
# =========================

geiger_results <- arrange(geiger_results, trait, aicc)

write_tsv(geiger_results,
          "model_fitting_geiger_results_archaea.txt")

# =========================
# SELECT BEST MODEL PER TRAIT
# =========================

best_models <- geiger_results %>%
  group_by(trait) %>%
  slice_min(aicc, n = 1) %>%
  ungroup()

get_param <- function(row) {
  model <- row$model
  if (model == "lambda") return(row$lambda)
  if (model == "kappa")  return(row$kappa)
  if (model == "delta")  return(row$delta)
  if (model == "EB")     return(row$a)
  return(NA_real_)
}

best_models <- best_models %>%
  rowwise() %>%
  mutate(param = get_param(cur_data())) %>%
  ungroup()

write_tsv(best_models,
          "best_model_per_trait_archaea.txt")

# =========================
# ANCESTRAL STATE RECONSTRUCTION
# =========================

dir.create("ASR_plots", showWarnings = FALSE)

ancestral_states <- list()
rescaled_trees <- list()
contmaps <- list()

for (trait_name in names(tip_values)) {
  
  message("==Processing Trait:", trait_name, "==")
  
  trait_vector <- tip_values[[trait_name]]
  
  best_model <- best_models %>%
    filter(trait == trait_name)
  
  t_sub <- drop.tip(tree2,
                    setdiff(tree2$tip.label,
                            names(trait_vector)))
  
  # rescale tree
  if (best_model$model %in% c("lambda","kappa","delta")) {
    t_rescaled <- phytools::rescale(
      t_sub,
      model = best_model$model,
      best_model$param
    )
  } else if (best_model$model == "EB") {
    t_rescaled <- phytools::rescale(
      t_sub,
      model = "EB",
      best_model$param
    )
  } else {
    t_rescaled <- t_sub
  }
  
  rescaled_trees[[trait_name]] <- t_rescaled
  
  # ancestral reconstruction (correct function)
  anc <- fastAnc(t_rescaled,
                 trait_vector,
                 vars = TRUE,
                 CI = TRUE)
  
  ancestral_states[[trait_name]] <- anc
  
  # continuous trait mapping
  cmap <- phytools::contMap(t_rescaled,
                            trait_vector,
                            plot = FALSE)
  
  pdf(paste0("ASR_plots/", trait_name, "_contMap.pdf"),
      width = 10,
      height = 50)
  
  plot(cmap,
       fsize = 0.3,
       lwd = 3,
       outline = FALSE,
       legend = TRUE)
  
  dev.off()
  
  contmaps[[trait_name]] <- cmap
}

saveRDS(ancestral_states, "ASR_archaea_states.rds")
saveRDS(rescaled_trees, "ASR_archaea_rescaled_trees.rds")
saveRDS(contmaps, "ASR_archaea_contmaps.rds")

message("==Pipeline complete! ==")

