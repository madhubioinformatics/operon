# scRepertoire: Single-Cell TCR/BCR Immune Profiling
Workflow for processing and analyzing single-cell T-cell receptor (TCR) and
B-cell receptor (BCR) data using scRepertoire in R. Designed to integrate
directly with Seurat and SingleCellExperiment workflows.

## Overview
scRepertoire processes filtered contig outputs from 10x Genomics Cell Ranger VDJ
pipeline and links clonotype information to single-cell RNA-seq expression data.
Supports 10x, AIRR, BD, MiXCR, TRUST4, and WAT3R single-cell clonal formats.

**Citation**:
- v2: Yang et al. (2025). scRepertoire 2: Enhanced and efficient toolkit for
  single-cell immune profiling. PLoS Computational Biology.
  https://doi.org/10.1371/journal.pcbi.1012760
- v1: Borcherding et al. (2020). scRepertoire: An R-based toolkit for
  single-cell immune receptor analysis. F1000Research.
  https://doi.org/10.12688/f1000research.22139.2

## Environment
- R (>= 4.1)
- Required R packages: scRepertoire, Seurat, BiocManager, immApex, ggplot2, dplyr

## Installation
```r
# Install from Bioconductor (stable)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scRepertoire")

# OR install latest from GitHub (requires immApex dependency)
remotes::install_github(c("BorchLab/immApex", "BorchLab/scRepertoire@devel"))
```

## Input Data
- Primary input: `filtered_contig_annotations.csv` from 10x Cell Ranger VDJ pipeline
  - Located at: `./outs/filtered_contig_annotations.csv`
- One CSV file per sample
- Paired scRNA-seq Seurat object (optional but recommended for integration)

## Workflow Steps

### 1. Load Contig Data
```r
library(scRepertoire)
library(Seurat)

# Load filtered contig files for each sample into a list
S1 <- read.csv("sample1/outs/filtered_contig_annotations.csv")
S2 <- read.csv("sample2/outs/filtered_contig_annotations.csv")
S3 <- read.csv("sample3/outs/filtered_contig_annotations.csv")

contig_list <- list(S1, S2, S3)

# Use built-in example data to test
data("contig_list")  # 8 samples, 4 patients, peripheral blood + BAL
```

### 2. Combine Contigs into Clonotypes

#### For TCR data:
```r
combined.TCR <- combineTCR(
  contig_list,
  samples = c("P17B", "P17L", "P18B", "P18L",
              "P19B", "P19L", "P20B", "P20L")
)
```

#### For BCR data:
```r
combined.BCR <- combineBCR(
  contig_list,
  samples = c("sample1", "sample2", "sample3"),
  call.related.clones = TRUE  # groups B cells by CDR3 similarity
)
```

### 3. Basic Clonal Analysis

#### Clonotype frequency and abundance
```r
# Number of unique clonotypes per sample
clonalQuant(combined.TCR, clone.call = "strict", chain = "both")

# Abundance of clonotypes (scaled vs unscaled)
clonalAbundance(combined.TCR, clone.call = "gene", scale = TRUE)

# CDR3 length distribution
clonalLength(combined.TCR, clone.call = "aa", chain = "TRB")
```

#### Clonal space homeostasis
```r
# How much space is occupied by rare vs expanded clones
clonalHomeostasis(combined.TCR, clone.call = "gene")

# Proportion of cells belonging to top N clones
clonalProportion(combined.TCR, clone.call = "gene",
                 clonalSplit = c(1, 5, 10, 100, 1000))
```

### 4. Repertoire Overlap Between Samples
```r
# Pairwise overlap heatmap (Jaccard, Morisita, overlap coefficient)
clonalOverlap(combined.TCR,
              clone.call = "strict",
              method = "morisita")
```

### 5. V/J Gene Usage Visualization
```r
# V gene usage bar chart
vizGenes(combined.TCR,
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "bar",
         scale = TRUE)

# Paired V-J gene usage heatmap
vizGenes(combined.TCR,
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap")
```

### 6. Diversity Analysis
```r
# Shannon, inverse Simpson, Chao1 diversity metrics
clonalDiversity(combined.TCR,
                clone.call = "gene",
                group.by = "sample",
                n.boots = 100)
```

### 7. Combine with Seurat Object
```r
# Load your Seurat scRNA-seq object
seu <- readRDS("seurat_object.rds")

# Add clonotype information to Seurat metadata by cell barcode
seu <- combineExpression(
  combined.TCR,
  seu,
  clone.call = "gene",
  group.by = "sample",
  proportion = TRUE
)

# Visualize clonal expansion on UMAP
DimPlot(seu, group.by = "cloneSize") +
  scale_color_manual(values = colorblind_vector(5))

# Highlight specific clonotypes on UMAP
highlightClonotypes(seu,
                    clone.call = "aa",
                    sequence = c("CASSLYYGYTF;CASSLAYYGYTF"))
```

### 8. Clonotype Tracking Across Conditions
```r
# Track clonotypes shared between peripheral blood and tumor
clonalNetwork(combined.TCR,
              clone.call = "strict",
              filter.clones = NULL,
              filter.identity = NULL)
```

## Conventions
- Always use `filtered_contig_annotations.csv` (not `all_contig_annotations.csv`) for cleaner data
- Use `clone.call = "strict"` for most stringent clonotype calling (gene + CDR3 aa)
- Use `clone.call = "aa"` for CDR3 amino acid sequence-based calling
- Save all figures to `figures/` directory with `ggsave()`
- Set `set.seed(42)` before any bootstrapping or diversity calculations
- Check cell barcode format consistency between Seurat object and contig files before `combineExpression()`

## Common Pitfalls
- Barcode prefixes/suffixes must match exactly between Seurat object and contig files — check with `head(colnames(seu))` vs `head(combined.TCR[[1]]$barcode)`
- `combineTCR()` and `combineBCR()` require separate calls — do not mix TCR and BCR contigs
- For multi-sample Seurat objects, cell barcodes must include sample prefix (e.g., `sample1_ACGT...`)
- TRUST4-reconstructed contigs from bulk RNA-seq need reformatting before loading into scRepertoire
- BCR somatic hypermutation means clonotype definitions are less straightforward than TCR — use `call.related.clones = TRUE` in `combineBCR()`

## References
- scRepertoire documentation: https://www.borch.dev/uploads/screpertoire/
- GitHub: https://github.com/BorchLab/scRepertoire
- Bioconductor: https://www.bioconductor.org/packages/release/bioc/html/scRepertoire.html
- v2 Citation: https://doi.org/10.1371/journal.pcbi.1012760
