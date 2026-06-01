# Stingless Bee Gut Microbiome Pipeline

Code repository for:

> **Gain and loss of gut symbionts during stingless bee diversification is linked to host group size**  
> Villabona et al. — *bioRxiv*

---

## Overview

A three-script R Markdown pipeline for 16S amplicon analysis of the gut microbiome of stingless bees (*Meliponini*), integrating absolute abundance estimates from qPCR and phylogenetic comparative methods.

```
01_Decontam_and_Filtering.Rmd   →   02_Microbiome_Analysis.Rmd   →   03_qPCR_and_AbsAbundance.Rmd
```

Each script reads inputs from `data/` and writes processed outputs back to `data/processed/` and figures to `figures/raw/`. Scripts must be run in order.

---

## Repository structure

```
SB_Pipeline/
├── data/
│   ├── raw/              # Input files (OTU table, taxonomy, metadata, phylogeny, qPCR)
│   └── processed/        # Intermediate outputs written by each script
├── figures/
│   ├── raw/              # PDF figures produced by scripts 02 and 03
│   └── final/            # Polished figures for publication
└── scripts/              # R Markdown scripts (run in numbered order)
```

---

## Scripts

### 01 — Decontam and Filtering
- Merges sequencing metadata with collection metadata
- Detects and removes contaminants using the prevalence method (`decontam`, threshold = 0.5)
- Filters samples with >25% contaminant load
- **Outputs:** `phyloseqOBJ.rds`, `muestras_a_eliminar.rds`, `contamdf.prev05.rds`, `metadata_2025.txt`

### 02 — Microbiome Analysis
- Alpha diversity (Shannon, Richness) with PGLS against colony size
- Beta diversity: CLR-PCA and Bray-Curtis PCoA with PERMANOVA
- Composition barplot annotated by biogeographic region
- Symbiont relative abundance by bee genus (boxplots)
- Pagel's lambda test and ancestral state reconstruction for symbiont presence/absence
- **Outputs:** `alpha_table.rds`, `pruned_tree.rds`, `metadata_symbionts.rds`, `discrete_*.rds`, figures

### 03 — qPCR and Absolute Abundance
- Processes qPCR triplicates: outlier filtering, efficiency correction, 16S copy number estimation
- Calculates absolute symbiont abundance = relative abundance × total 16S copies per µL
- Barplots and boxplots of absolute abundance by bee genus
- PGLS of absolute abundance against colony size
- **Outputs:** `qpcr_filtered.rds`, `unique_all_avg.rds`, `estimated_abs_abundance_*.tsv`, `OTUs_AbsAbun.tsv`, figures

---

## Dependencies

R ≥ 4.2. Install required packages:

```r
install.packages(c("tidyverse", "vegan", "ape", "nlme", "cowplot",
                   "RColorBrewer", "ggpubr", "gridExtra", "reshape2",
                   "vroom", "phytools", "corHMM", "smplot2", "funrar"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "decontam"))

# microViz (CRAN)
install.packages("microViz", repos = c(davidbarnett = "https://david-barnett.r-universe.dev",
                                        getOption("repos")))
```

---

## Running the pipeline

Set the working directory to the repository root before knitting each script:

```r
# In RStudio: Session > Set Working Directory > To Source File Location
# or:
setwd("path/to/SB_Pipeline")
rmarkdown::render("scripts/01_Decontam_and_Filtering.Rmd")
rmarkdown::render("scripts/02_Microbiome_Analysis.Rmd")
rmarkdown::render("scripts/03_qPCR_and_AbsAbundance.Rmd")
```

---

## Data

Raw input files are in `data/raw/`:

| File | Description |
|------|-------------|
| `OTU_table_nc.tsv` | DADA2 ASV table (16S V4, negative controls removed) |
| `Orb_Tax_table.txt` | Taxonomy table with manual LCA corrections |
| `metadata_table.txt` | Sequencing run metadata |
| `MappingFile_Meliponini_Claus.txt` | Collection metadata (Claus et al.) |
| `meliponini5genes.nex` | Five-gene Meliponini phylogeny (NEXUS) |
| `qPCR_All_triplicates.csv` | Raw qPCR triplicates for all samples and symbionts |
| `qPCR_All_Triplicates_WB.csv` | qPCR triplicates (whole-body extractions) |

---

## Contact

Nickole Villabona — nvillabo@uci.edu — University of California, Irvine
