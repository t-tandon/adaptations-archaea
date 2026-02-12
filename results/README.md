# Results Summary

## 1. Exploratory Patterns

Across major archaeal phyla:

* **frac_ilvwygerkp** shows a strong positive relationship with OGT in nearly all clades.
* **ratio_charged_polar** generally increases with OGT.
* **tRNA GC content** shows a positive association in many thermophilic groups.
* Whole-genome GC content and genome size show weaker trends.

Class-level analysis shows:

* Strong clustering by class.
* More modest/noisy patterns in small classes.
* **frac_ilvwygerkp** still gives consistent positive slopes


## 2. Phylogenetic Signal in Traits


| Trait               | Best Model | Key Parameter | Interpretation                             |
| ------------------- | ---------- | ------------- | ------------------------------------------ |
| OGT                 | λ          | λ ≈ 0.99      | Very strong phylogenetic signal       |
| Genome size         | λ          | λ ≈ 0.88      | Strong phylogenetic signal             |
| frac_ilvwygerkp     | δ          | δ ≈ 0.19      | Evolution in deeper nodes |
| GC content          | δ          | δ ≈ 0.71      | Somewhat time-dependent evolution        |
| ratio_charged_polar | EB         | a ≈ -0.92     | Early burst–like pattern                   |
| trna_gc             | δ          | δ ≈ 0.028     | Early divergence                    |


* **OGT is highly phylogenetically conserved** (λ ≈ 0.99).
* Many molecular traits show early divergence patterns, which might mean rapid early adaptation followed by constraint.


## 3. Per-Phylum GLS Results

Across phyla with enough statistical power:

### Universally Strong Predictor: **frac_ilvwygerkp**

* Positive and highly significant in nearly every major phylum.
* Represents hydrophobic and thermophily-enriched amino acids.

### Secondary Predictors

#### **ratio_charged_polar**

* Frequently positive and significant.
* Supports electrostatic stabilization mechanisms.

#### **tRNA GC content**

* Significant in multiple thermophilic clades.
* Likely contributes to RNA structural stability at high temperature.


### Weak or Inconsistent Predictors

* Genome size
* Whole-genome GC content

These variables sometimes show clade-specific and are not consistent predictors of thermophily.


## 4. Per-Class Analysis

At the class level:

* ILVWYGERKP is still a strong predictor in thermophilic classes.
* Even with lower power, smaller classes show positive slopes
Indicates that proteome composition is the main predictor of OGT. 




