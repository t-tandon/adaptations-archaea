# adaptations-archaea

## project workflow
1. Select high-quality representative genomes
   - Completeness: >90%
   - Contamination: <5%
   - One per genus
3. Build phylogenetic tree
  - Identify universal marker proteins with Prokka
  - Align, trim, concatenate
  - Infer tree with IQ-Tree
3. Calculate genomic features for each genome
  - Genome size, gene density, GC3, intergenic lengths, amino acid composition, etc.
4. Test trait correlation with OGT in R
  - Using phylogenetic generalized least square (PGLS) or Phylogenetic Independent Contrasts (PIC) to control for shared ancestry
5. Visualize features across tree to see lineage-specific patterns. 
