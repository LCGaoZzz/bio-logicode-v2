# Tool-Specific Patterns Reference

## Overview

This reference documents common code patterns for bioinformatics tools, helping Bio-LogiCode recognize and map tool usage to biological intent.

## R / Seurat

### Loading Data

**scRNA-seq:**
```r
library(Seurat)

# Standard 10x data
counts <- Read10X("filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts, project = "project_name")

# H5 format
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
seurat_obj <- CreateSeuratObject(counts)
```

**Visium spatial:**
```r
# Load Visium data with spatial coordinates
visium_obj <- Load10X_Spatial(
  data.dir = "outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial"
)
```

**Xenium spatial:**
```r
# Seurat v5+
xenium_obj <- LoadXenium("path/to/xenium_output/", fov = "fov")
```

### QC and Filtering

```r
# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Filtering
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                              nFeature_RNA < 5000 &
                              nCount_RNA > 500 &
                              percent.mt < 20)
```

**Key patterns to recognize:**
- `PercentageFeatureSet()` - Calculate QC metrics
- `subset()` - Filter cells/spots
- `nFeature_RNA`, `nCount_RNA`, `percent.mt` - Common QC metrics

### Normalization

```r
# Log-normalization (standard)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# SCTransform (variance-stabilizing)
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt")
```

**Key patterns:**
- `NormalizeData()` - Log-normalization
- `SCTransform()` - Variance-stabilizing transformation (recommended for spatial)

### Dimensionality Reduction

```r
# PCA
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# tSNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
```

**Key patterns:**
- `FindVariableFeatures()` - Select highly variable genes
- `ScaleData()` - Z-score scaling
- `RunPCA()`, `RunUMAP()`, `RunTSNE()` - Dimensionality reduction

### Clustering

```r
# Graph-based clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
```

**Key parameters to extract:**
- `dims` - Which PCs used
- `resolution` - Clustering granularity (higher = more clusters)

### Differential Expression

```r
# Find markers for all clusters
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Find markers for specific cluster
cluster1_markers <- FindMarkers(seurat_obj, ident.1 = 1, min.pct = 0.25)

# Compare two groups
de_genes <- FindMarkers(seurat_obj, ident.1 = "group1", ident.2 = "group2")
```

**Key parameters:**
- `only.pos` - Only positive markers
- `min.pct` - Minimum fraction of cells expressing
- `logfc.threshold` - Minimum log fold-change

### Spatial Analysis (Seurat)

```r
# Spatially variable features
seurat_obj <- FindSpatiallyVariableFeatures(seurat_obj,
                                            assay = "Spatial",
                                            features = VariableFeatures(seurat_obj),
                                            selection.method = "markvariogram")

# Spatial visualization
SpatialFeaturePlot(seurat_obj, features = "gene_name")
SpatialDimPlot(seurat_obj, label = TRUE)
```

**Key patterns:**
- `FindSpatiallyVariableFeatures()` - Identify spatially variable genes
- `SpatialFeaturePlot()`, `SpatialDimPlot()` - Spatial visualization

## Python / Scanpy

### Loading Data

```python
import scanpy as sc

# Standard 10x data
adata = sc.read_10x_mtx("filtered_feature_bc_matrix/")

# H5 format
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")

# H5AD format (AnnData)
adata = sc.read_h5ad("data.h5ad")

# Visium spatial
adata = sc.read_visium("outs/")

# Loom format
adata = sc.read_loom("data.loom")
```

### QC and Filtering

```python
# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# OR manually
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
```

**Key patterns:**
- `sc.pp.calculate_qc_metrics()` - Calculate QC metrics
- `sc.pp.filter_cells()`, `sc.pp.filter_genes()` - Filter cells/genes
- `pct_counts_mt` - Mitochondrial percentage

### Normalization

```python
# Log-normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
```

**Key patterns:**
- `sc.pp.normalize_total()` - Total-count normalization
- `sc.pp.log1p()` - Log transformation
- `sc.pp.highly_variable_genes()` - Select highly variable genes

### Dimensionality Reduction

```python
# PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# UMAP
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# tSNE
sc.tl.tsne(adata, n_pcs=40)
```

**Key patterns:**
- `sc.pp.scale()` - Z-score scaling
- `sc.tl.pca()` - PCA
- `sc.pp.neighbors()` - Build k-nearest neighbor graph
- `sc.tl.umap()`, `sc.tl.tsne()` - Dimensionality reduction

### Clustering

```python
# Leiden clustering (recommended)
sc.tl.leiden(adata, resolution=0.8)

# Louvain clustering
sc.tl.louvain(adata, resolution=0.8)
```

**Key parameters:**
- `resolution` - Clustering granularity

### Differential Expression

```python
# Rank genes for all clusters
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Methods: 'wilcoxon', 't-test', 'logreg'
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

# Get results
result = adata.uns['rank_genes_groups']
```

**Key parameters:**
- `method` - Statistical test ('wilcoxon', 't-test', 'logreg')

## Python / Squidpy (Spatial Analysis)

### Loading Spatial Data

```python
import squidpy as sq

# Visium
adata = sq.read.visium("outs/")

# Xenium
adata = sq.read.xenium("path/to/xenium_output/")
```

### Spatial Neighbor Graph

```python
# Build spatial neighbor graph
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

# For Visium (grid-based)
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
```

**Key parameters:**
- `coord_type` - "generic" (Xenium) or "grid" (Visium)
- `n_neighs` - Number of neighbors

### Spatial Autocorrelation

```python
# Moran's I (spatial autocorrelation)
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=adata.var_names[:100],
    n_perms=100,
    n_jobs=1
)

# Geary's C
sq.gr.spatial_autocorr(adata, mode="geary")
```

**Key patterns:**
- `sq.gr.spatial_autocorr()` - Moran's I, Geary's C
- `mode="moran"` or `mode="geary"`

### Neighborhood Enrichment

```python
# Cell type co-localization
sq.gr.nhood_enrichment(adata, cluster_key="leiden")

# Visualization
sq.pl.nhood_enrichment(adata, cluster_key="leiden")
```

**Key patterns:**
- `sq.gr.nhood_enrichment()` - Test for cell type co-localization
- `cluster_key` - Column in `adata.obs` with cell type labels

### Ligand-Receptor Interaction Analysis

```python
# Ligand-receptor analysis
sq.gr.ligrec(
    adata,
    n_perms=100,
    cluster_key="leiden",
    copy=False,
    use_raw=False
)

# Visualization
sq.pl.ligrec(adata, cluster_key="leiden")
```

**Key patterns:**
- `sq.gr.ligrec()` - Ligand-receptor interaction analysis
- `n_perms` - Number of permutations for significance testing

### Spatially Variable Genes

```python
# SpatialDE (external package often used with Squidpy)
import SpatialDE

# Prepare data
counts = adata.to_df()
sample_info = adata.obs[['x', 'y']]

# Run SpatialDE
results = SpatialDE.run(counts, sample_info)
```

## Python / SpatialData

### Loading Xenium Data

```python
import spatialdata as sd

# Load Xenium data (comprehensive spatial data structure)
sdata = sd.read_xenium("path/to/xenium_output/")

# Access components
transcripts = sdata.points['transcripts']  # Transcript coordinates
cells = sdata.shapes['cells']  # Cell boundaries
images = sdata.images['morphology']  # Morphology images
table = sdata.table  # AnnData table (cell × gene)
```

**Key patterns:**
- `sd.read_xenium()` - Load Xenium data into SpatialData structure
- `.points`, `.shapes`, `.images`, `.table` - Access different data modalities

## R / DESeq2 (Bulk RNA-seq)

### Standard Workflow

```r
library(DESeq2)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ condition)

# Pre-filtering (optional)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "treated", "control"))

# Shrink log2 fold changes
res <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")
```

**Key parameters to extract:**
- `design` - Model formula (e.g., `~ condition`, `~ batch + condition`)
- `contrast` - Which comparison (e.g., treated vs control)
- `lfcShrink` type - "apeglm", "ashr", "normal"

### Filtering Results

```r
# Filter by adjusted p-value and log2 fold change
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
```

**Key thresholds:**
- `padj` - Adjusted p-value (FDR)
- `log2FoldChange` - Effect size

## R / edgeR (Bulk RNA-seq)

### Standard Workflow

```r
library(edgeR)

# Create DGEList object
dge <- DGEList(counts = counts, group = group)

# Filtering
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalization
dge <- calcNormFactors(dge)

# Design matrix
design <- model.matrix(~ group)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit model
fit <- glmQLFit(dge, design)

# Test
qlf <- glmQLFTest(fit, coef = 2)
res <- topTags(qlf, n = Inf)
```

**Key parameters:**
- `design` - Model matrix
- `coef` - Which coefficient to test
- FDR threshold in `topTags()`

## Common Pattern Recognition Rules

### Batch Correction

**R / Seurat:**
```r
# Harmony
library(harmony)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch")

# Seurat integration
seurat_list <- SplitObject(seurat_obj, split.by = "batch")
anchors <- FindIntegrationAnchors(object.list = seurat_list)
seurat_integrated <- IntegrateData(anchorset = anchors)
```

**Python / Scanpy:**
```r
# Combat
sc.pp.combat(adata, key='batch')

# Harmony
import harmonypy
ho = harmonypy.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
adata.obsm['X_pca_harmony'] = ho.Z_corr.T

# scVI
import scvi
scvi.model.SCVI.setup_anndata(adata, batch_key='batch')
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scVI'] = model.get_latent_representation()
```

**Key patterns to recognize:**
- `RunHarmony()`, `harmonypy.run_harmony()` - Harmony batch correction
- `sc.pp.combat()` - ComBat batch correction
- `scvi.model.SCVI` - scVI deep learning batch correction

### Cell Type Annotation

**R / Seurat:**
```r
# Manual annotation
new.cluster.ids <- c("CD4 T", "CD8 T", "B cell", "Monocyte", ...)
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# SingleR (reference-based)
library(SingleR)
pred <- SingleR(test = seurat_obj, ref = ref_data, labels = ref_data$label)
```

**Python / Scanpy:**
```python
# Manual annotation
cluster_to_celltype = {
    '0': 'CD4 T cell',
    '1': 'CD8 T cell',
    '2': 'B cell',
    ...
}
adata.obs['celltype'] = adata.obs['leiden'].map(cluster_to_celltype)

# Automated annotation (celltypist)
import celltypist
model = celltypist.models.Model.load(model='Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model=model)
```

**Key patterns:**
- Manual mapping: cluster IDs → cell type labels
- Reference-based: SingleR, celltypist, scmap

## Detection Strategy Summary

When analyzing code, look for these import/library patterns:

**R:**
- `library(Seurat)` → scRNA-seq or spatial (Visium/Xenium)
- `library(DESeq2)` → Bulk RNA-seq DE
- `library(edgeR)` → Bulk RNA-seq DE
- `library(harmony)` → Batch correction
- `library(SingleR)` → Cell type annotation

**Python:**
- `import scanpy as sc` → scRNA-seq
- `import squidpy as sq` → Spatial analysis
- `import spatialdata as sd` → Xenium spatial data
- `import scvi` → Deep learning (batch correction, integration)
- `import SpatialDE` → Spatially variable genes

Then trace function calls to map biological intent → computational method → code evidence.
