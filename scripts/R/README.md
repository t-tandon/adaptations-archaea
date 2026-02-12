# README

This folder contains the R scripts used to perform phylogenetically informed regression analyses on prokaryotic genomic traits.

The analysis was conducted **twice at different taxonomic resolutions**:

* **Genus-level representatives (~700 sequences)**
* **Species-level representatives (~3000 sequences)**

Initial analyses attempted to use **Phylogenetic Generalized Least Squares (PGLS)**.
However, when working with the **larger species-level dataset**, PGLS frequently failed due to convergence issues in estimating phylogenetic signal (λ)

To ensure robust and comparable analyses across taxonomic levels, the workflow was adapted to use **Generalized Least Squares (GLS) with Pagel’s λ correlation structure**

---


