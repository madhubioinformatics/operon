# SeuratExtend: Enhanced Toolkit for scRNA-seq Analysis
Extended workflow for advanced scRNA-seq analysis and visualization using SeuratExtend,
built upon Seurat. Covers enhanced visualization, GSEA, trajectory analysis, regulatory
networks, and Python tool integration within R.

## Overview
SeuratExtend expands Seurat with advanced visualization, integrated pathway analysis,
seamless Python tool integration (scVelo, SCENIC, Palantir), and utility functions.
Designed to produce publication-ready figures while remaining beginner-friendly.

**Citation**: Hua et al. (2025). SeuratExtend: streamlining single-cell RNA-seq analysis
through an integrated and intuitive framework. *GigaScience* 14, giaf076.
https://doi.org/10.1093/gigascience/giaf076

## Environment
- R (>= 4.1)
- Required R packages: Seurat, SeuratExtend, ggplot2
- Optional Python env for trajectory tools: scVelo, Palantir, MAGIC, CellRank, pySCENIC

## Installation
```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("huayc09/SeuratExtend")
```

## Workflow Steps

### 1. Standard Seurat Pipeline (Prerequisite)
Run the standard Seurat workflow before using SeuratExtend:
```r
library(Seurat)
library(SeuratExtend)
```

### 2. Enhanced Visualization

#### UMAP with Arrow Theme
```r
DimPlot2(seu, theme = theme_umap_arrows())
```

#### Cluster Distribution Across Samples
```r
ClusterDistrBar(seu$orig.ident, seu$cluster)
```

#### Heatmap of Marker Genes (z-score)
```r
genes.zscore <- CalcStats(
  seu,
  features = VariableFeatures(seu),
  group.by = "cluster",
  order = "p",
  n = 4
)
Heatmap(genes.zscore, lab_fill = "zscore")
```

#### Enhanced Dot Plot with Grouped Markers
```r
grouped_features <- list(
  "B_cell_markers"  = c("MS4A1", "CD79A"),
  "T_cell_markers"  = c("CD3D", "CD8A", "IL7R"),
  "Myeloid_markers" = c("CD14", "FCGR3A", "S100A8")
)
DotPlot2(seu, features = grouped_features)
```

#### Violin Plot with Statistics
```r
VlnPlot2(
  seu,
  features = c("CD3D", "CD14", "CD79A"),
  group.by = "cluster",
  cells = WhichCells(seu, idents = c("B cell", "CD8 T cell", "Mono CD14")),
  stat.method = "wilcox.test"
)
```

#### Multi-Marker UMAP (3 markers, RYB coloring)
```r
FeaturePlot3(seu, feature.1 = "CD3D", feature.2 = "CD14", feature.3 = "CD79A", pt.size = 1)
```

#### Volcano Plot
```r
VolcanoPlot(seu, ident.1 = "B cell", ident.2 = "CD8 T cell")
```

### 3. Gene Set Enrichment Analysis (GSEA)

#### GO or Reactome Database
```r
options(spe = "human")  # or "mouse"
seu <- GeneSetAnalysisGO(seu, parent = "immune_system_process", n.min = 5)
matr <- RenameGO(seu@misc$AUCell$GO$immune_system_process)
go_zscore <- CalcStats(matr, f = seu$cluster, order = "p", n = 3)
Heatmap(go_zscore, lab_fill = "zscore")
```

#### GSEA Plot for a Specific Pathway
```r
GSEAplot(
  seu,
  ident.1 = "B cell",
  ident.2 = "CD8 T cell",
  title = "GO:0042113 B cell activation",
  geneset = GO_Data$human$GO2Gene[["GO:0042113"]]
)
```

#### Search and Filter Pathways
```r
SearchPathways(term = "interferon", database = "GO")
FilterGOTerms(parent = "immune_system_process", n.min = 5, n.max = 500)
```

### 4. Trajectory and Pseudotime Analysis

#### Palantir (Pseudotime in R)
```r
# Compute diffusion map
seu <- Palantir.RunDM(seu)
DimPlot2(seu, reduction = "ms")

# Compute pseudotime (set start cell manually)
seu <- Palantir.Pseudotime(seu, start_cell = "your_start_cell_barcode")

# Visualize gene trends along trajectory
ps <- seu@misc$Palantir$Pseudotime
GeneTrendCurve.Palantir(seu, pseudotime.data = ps, features = c("CD14", "FCGR3A"))
GeneTrendHeatmap.Palantir(seu, features = VariableFeatures(seu)[1:20],
                          pseudotime.data = ps, lineage = "fate1")
```

#### scVelo (RNA Velocity via Python)
```r
scVelo.SeuratToAnndata(
  seu,
  filename = "mye_small.h5ad",
  velocyto.loompath = "pbmc10k_mye_small.loom",
  prefix = "sample1_",
  postfix = "-1"
)
scVelo.Plot(color = "cluster", basis = "ms_cell_embeddings",
            save = "scvelo_output.png", figsize = c(5, 4))
```

#### MAGIC (Denoising Gene Expression)
```r
seu <- Palantir.Magic(seu)
```

### 5. Gene Regulatory Networks (SCENIC)

#### Import pySCENIC Results into Seurat
```r
seu <- ImportPyscenicLoom("pyscenic_integrated-output.loom", seu = seu)

# Visualize regulon activity alongside gene expression
DimPlot2(
  seu,
  features = c("cluster", "CEBPA", "tf_CEBPA"),
  cols = list("tf_CEBPA" = "OrRd"),
  theme = NoAxes()
)

# Waterfall plot comparing regulon activity between cell types
DefaultAssay(seu) <- "TF"
WaterfallPlot(seu, features = rownames(seu),
              ident.1 = "Mono CD14", ident.2 = "CD8 T cell",
              exp.transform = FALSE, top.n = 20)
```

### 6. Utility Functions

#### Gene Name Conversions
```r
HumanToMouseGenesymbol(c("CD3D", "CD14", "MS4A1"))
MouseToHumanGenesymbol(c("Cd3d", "Cd14"))
EnsemblToGenesymbol(c("ENSG00000167286"), spe = "human")
```

#### Cluster Statistics
```r
CalcStats(seu, features = c("CD3D", "CD14"), group.by = "cluster")
feature_percent(seu, features = c("CD3D", "CD14"), group.by = "cluster")
```

#### Run Full Standard Seurat Pipeline in One Step
```r
seu <- RunBasicSeurat(seu)
```

## Conventions
- Always run standard Seurat pipeline first before SeuratExtend functions
- Set up Python conda environment once: `create_condaenv_seuratextend()` for trajectory tools
- Save plots to `figures/` folder with consistent naming
- Use `random_state = 42` for reproducibility in Python-integrated tools
- For multi-sample datasets, check batch effects with `ClusterDistrBar()` before downstream analysis

## Common Pitfalls
- Python tools (scVelo, Palantir, SCENIC) require a properly configured conda environment — run `create_condaenv_seuratextend()` before first use
- `options(spe = "human")` or `"mouse"` must be set before GSEA functions
- SCENIC loom file must be generated externally with pySCENIC before importing with `ImportPyscenicLoom()`
- For `Palantir.Pseudotime()`, the start cell barcode must exactly match a cell name in your Seurat object
- `FeaturePlot3` displays exactly 3 features using RYB blending — use `FeaturePlot3.grid` for more

## References
- SeuratExtend GitHub: https://github.com/huayc09/SeuratExtend
- Full tutorial: https://huayc09.github.io/SeuratExtend/index.html
- Hua et al. (2025) GigaScience: https://doi.org/10.1093/gigascience/giaf076
