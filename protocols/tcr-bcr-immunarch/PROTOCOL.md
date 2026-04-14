# immunarch: TCR/BCR Immune Repertoire Analysis
Comprehensive workflow for analyzing T-cell receptor (TCR) and B-cell receptor (BCR)
repertoire data using immunarch in R. Supports bulk and single-cell formats including
10x Genomics, AIRR, MiXCR, IMGT, ImmunoSEQ, and more.

## Overview
immunarch is an R package for painless analysis of adaptive immune receptor repertoire
(AIRR) data. It automatically detects file formats, supports all popular TCR/BCR
sequencing platforms, and is fully interoperable with Seurat for single-cell workflows.

**Citation**: ImmunoMind Team. (2019). immunarch: An R Package for Painless
Bioinformatics Analysis of T-Cell and B-Cell Immune Repertoires. Zenodo.
http://doi.org/10.5281/zenodo.3367200

## Environment
- R (>= 3.5)
- Required R packages: immunarch, pak, ggplot2, dplyr

## Installation
```r
# Install pak first (recommended installer)
install.packages("pak", repos = sprintf(
  "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
  .Platform$pkgType, R.Version()$os, R.Version()$arch))

# Install latest stable version
pak::pkg_install("immunomind/immunarch@0.9.1")

# OR install from CRAN
install.packages("immunarch")
```

## Input Data Formats Supported
- 10x Genomics (scTCR/BCR)
- AIRR-seq
- MiXCR, MiGEC, MiGMap
- IMGT, ImmunoSEQ
- VDJtools, ArcherDX
- R data frames and data.table

## Workflow Steps

### 1. Load Data
```r
library(immunarch)

# Load test dataset to explore immunarch
data(immdata)

# Load your own data (auto-detects format)
immdata <- repLoad("path/to/your/folder/with/repertoires")

# immdata is a list with two elements:
# immdata$data  — list of repertoire data frames (one per sample)
# immdata$meta  — metadata table with sample information
```

### 2. Explore Basic Statistics
```r
# CDR3 length distribution
repExplore(immdata$data, "lens") %>% vis()

# Number of clonotypes per sample
repExplore(immdata$data, "counts") %>% vis()

# Gene usage (V, J, D segments)
repExplore(immdata$data, "volume") %>% vis()
```

### 3. Clonality Analysis
```r
# Clonal space homeostasis (relative abundance of clonotypes)
repClonality(immdata$data, "homeo") %>% vis()

# Top clones proportion
repClonality(immdata$data, "top", .head = c(10, 100, 1000)) %>% vis()

# Rare clones analysis
repClonality(immdata$data, "rare") %>% vis()
```

### 4. Repertoire Overlap Between Samples
```r
# Compute pairwise overlap
ov <- repOverlap(immdata$data)
vis(ov)

# Cluster samples by overlap similarity
ov.kmeans <- repOverlapAnalysis(ov, .method = "mds+kmeans")
vis(ov.kmeans)
```

### 5. V/J Gene Usage Analysis
```r
# V gene usage for all samples
gu <- geneUsage(immdata$data)
vis(gu)

# Compare V gene usage by clinical status
gu_grouped <- geneUsage(immdata$data)
vis(gu_grouped, .by = "Status", .meta = immdata$meta)

# Cluster samples by V gene usage similarity
gu.clust <- geneUsageAnalysis(gu, .method = "js+hclust")
vis(gu.clust)
```

### 6. Repertoire Diversity
```r
# Shannon entropy diversity
repDiversity(immdata$data, "shannon") %>%
  vis(.by = "Status", .meta = immdata$meta)

# Chao1 diversity estimate
repDiversity(immdata$data, "chao1") %>%
  vis(.by = "Status", .meta = immdata$meta)

# D50 index
repDiversity(immdata$data, "d50") %>% vis()
```

### 7. Clonotype Tracking Across Samples
```r
# Track specific clonotypes across timepoints or conditions
tc <- trackClonotypes(immdata$data, list(1, "aa"), .col = "aa")
vis(tc)
```

### 8. K-mer and Sequence Analysis
```r
# K-mer frequency analysis of CDR3 sequences
km <- getKmers(immdata$data[[1]], 3)
vis(km)
```

### 9. Integration with Seurat (Single-Cell)
```r
# After running standard Seurat workflow on mRNA data:
library(Seurat)
library(immunarch)

# Load TCR contigs from 10x Cell Ranger VDJ output
immdata <- repLoad("path/to/vdj/output/")

# Merge clonotype data into Seurat metadata by cell barcode
# Match barcodes between Seurat object and immunarch data
```

## Conventions
- Always inspect `immdata$meta` to confirm sample metadata is correctly loaded
- Use `vis()` after every analysis function for quick visualization
- Save results with `ggsave("figures/plot_name.pdf", width=8, height=6)`
- For multi-sample comparisons, always pass `.meta = immdata$meta` and `.by = "Group"` to group by condition
- Use `select()` from dplyr to subset samples before analysis

## Common Pitfalls
- File format is auto-detected but folder must contain only repertoire files — remove unrelated files
- For 10x single-cell data, ensure barcodes match exactly between VDJ and GEX libraries
- Clonality metrics assume sufficient sequencing depth — check `repExplore(immdata$data, "counts")` first
- B-cell clonotype definitions differ from T-cells due to somatic hypermutation — interpret BCR clonality with caution
- Always set `random_state` or `set.seed(42)` before any clustering steps

## References
- immunarch website: https://immunarch.com
- GitHub: https://github.com/immunomind/immunarch
- Citation: http://doi.org/10.5281/zenodo.3367200
