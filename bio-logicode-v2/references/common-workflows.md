# Common Bioinformatics Workflows Reference

## Overview

This reference documents standard workflow patterns for common bioinformatics analyses, helping Bio-LogiCode recognize typical analysis pipelines and their expected steps.

## scRNA-seq (Single-Cell RNA-seq) Workflow

### Standard Pipeline Steps

1. **Quality Control (QC)**
   - **Biological intent:** Remove low-quality cells and doublets
   - **Common metrics:**
     - Min genes per cell (typical: 200-500)
     - Max genes per cell (typical: 2500-5000, removes doublets)
     - Min UMI per cell (typical: 500-1000)
     - Max mitochondrial percentage (typical: 5-20%)
   - **Tools:** Seurat `subset()`, Scanpy `sc.pp.filter_cells()`

2. **Normalization**
   - **Biological intent:** Account for sequencing depth differences between cells
   - **Common methods:**
     - Log-normalization (scale to 10,000 UMI per cell, then log1p)
     - SCTransform (variance-stabilizing transformation)
     - scran pooling-based normalization
   - **Tools:** Seurat `NormalizeData()` or `SCTransform()`, Scanpy `sc.pp.normalize_total()` + `sc.pp.log1p()`

3. **Feature Selection**
   - **Biological intent:** Identify highly variable genes for downstream analysis
   - **Common parameters:**
     - Number of features: 2000-3000
     - Selection method: vst, mean.var.plot, dispersion
   - **Tools:** Seurat `FindVariableFeatures()`, Scanpy `sc.pp.highly_variable_genes()`

4. **Dimensionality Reduction**
   - **Biological intent:** Reduce data to principal components capturing biological variation
   - **Common steps:**
     - Scaling (z-score normalization)
     - PCA (typically 30-50 PCs)
     - UMAP or tSNE for visualization
   - **Tools:** Seurat `ScaleData()` + `RunPCA()` + `RunUMAP()`, Scanpy `sc.pp.scale()` + `sc.tl.pca()` + `sc.tl.umap()`

5. **Batch Correction (if multiple batches)**
   - **Biological intent:** Remove technical batch effects while preserving biological variation
   - **Common methods:**
     - Harmony
     - ComBat
     - Seurat integration
     - scVI
   - **Tools:** Harmony `RunHarmony()`, Scanpy `sc.pp.combat()`, scVI

6. **Clustering**
   - **Biological intent:** Group cells with similar expression profiles
   - **Common methods:**
     - Graph-based clustering (Louvain, Leiden)
     - K-means
   - **Common parameters:**
     - Resolution: 0.4-1.2 (higher = more clusters)
     - Number of PCs: 20-50
   - **Tools:** Seurat `FindNeighbors()` + `FindClusters()`, Scanpy `sc.pp.neighbors()` + `sc.tl.leiden()`

7. **Differential Expression (DE)**
   - **Biological intent:** Identify marker genes for each cluster or compare conditions
   - **Common methods:**
     - Wilcoxon rank-sum test
     - t-test
     - MAST (for scRNA-seq)
     - DESeq2 (pseudobulk)
   - **Common thresholds:**
     - Adjusted p-value (FDR): < 0.05
     - Log fold-change: > 0.25 or > 0.5
     - Min percentage expressing: > 0.25
   - **Tools:** Seurat `FindMarkers()` / `FindAllMarkers()`, Scanpy `sc.tl.rank_genes_groups()`

8. **Cell Type Annotation**
   - **Biological intent:** Assign biological cell type labels to clusters
   - **Common approaches:**
     - Manual annotation (marker genes)
     - Reference-based (SingleR, celltypist)
     - Automated (scmap, Garnett)
   - **Tools:** Seurat `RenameIdents()`, SingleR, celltypist

9. **Trajectory Analysis (optional)**
   - **Biological intent:** Infer developmental trajectories or differentiation paths
   - **Common tools:**
     - Monocle3
     - Slingshot
     - PAGA
     - RNA velocity
   - **Output:** Pseudotime ordering, lineage assignments

10. **Visualization**
    - **Common plots:**
      - UMAP/tSNE with clusters
      - Feature plots (gene expression)
      - Violin plots (marker genes)
      - Heatmaps (top markers)
      - Dot plots (marker expression across clusters)

### Typical Parameter Ranges

| Step | Parameter | Typical Range | Notes |
|------|-----------|---------------|-------|
| QC | min_genes | 200-500 | Lower for sparse data |
| QC | max_genes | 2500-5000 | Removes doublets |
| QC | max_mt% | 5-20% | Higher for certain tissues |
| Normalization | scale_factor | 10000 | Standard for log-norm |
| Feature selection | n_features | 2000-3000 | More for complex datasets |
| PCA | n_pcs | 30-50 | Check elbow plot |
| Clustering | resolution | 0.4-1.2 | Higher = more clusters |
| DE | min.pct | 0.1-0.25 | Min fraction expressing |
| DE | logfc.threshold | 0.25-0.5 | Min log fold-change |

## Bulk RNA-seq Workflow

### Standard Pipeline Steps

1. **Quality Control (Raw Reads)**
   - **Biological intent:** Assess sequencing quality
   - **Tools:** FastQC, MultiQC
   - **Metrics:** Per-base quality, adapter content, duplication rates

2. **Adapter Trimming (if needed)**
   - **Biological intent:** Remove adapter sequences
   - **Tools:** Cutadapt, Trim Galore, fastp

3. **Alignment**
   - **Biological intent:** Map reads to reference genome
   - **Common tools:**
     - STAR (splice-aware, fast)
     - HISAT2 (splice-aware, memory-efficient)
     - Salmon/Kallisto (pseudo-alignment, fast)
   - **Output:** BAM files or transcript counts

4. **Quantification**
   - **Biological intent:** Count reads per gene/transcript
   - **Common tools:**
     - featureCounts (from BAM)
     - RSEM (from BAM)
     - Salmon/Kallisto (direct quantification)
   - **Output:** Count matrix (genes × samples)

5. **Quality Control (Count Matrix)**
   - **Biological intent:** Assess sample quality and detect outliers
   - **Common checks:**
     - Total read counts per sample
     - Number of detected genes
     - PCA to detect batch effects or outliers
   - **Tools:** DESeq2 `plotPCA()`, edgeR `plotMDS()`

6. **Filtering**
   - **Biological intent:** Remove lowly expressed genes
   - **Common thresholds:**
     - Min counts: 10-20 across all samples
     - Min CPM: 1-2 in at least N samples (N = smallest group size)
   - **Tools:** DESeq2 `rowSums()`, edgeR `filterByExpr()`

7. **Normalization**
   - **Biological intent:** Account for library size and composition differences
   - **Common methods:**
     - DESeq2: Median of ratios
     - edgeR: TMM (Trimmed Mean of M-values)
     - limma-voom: TMM + voom transformation
   - **Tools:** DESeq2 `DESeq()`, edgeR `calcNormFactors()`

8. **Differential Expression**
   - **Biological intent:** Identify genes differentially expressed between conditions
   - **Common methods:**
     - DESeq2 (negative binomial)
     - edgeR (negative binomial)
     - limma-voom (linear modeling)
   - **Common thresholds:**
     - Adjusted p-value (FDR): < 0.05 or < 0.01
     - Log2 fold-change: > 1 (2-fold change)
   - **Tools:** DESeq2 `results()`, edgeR `glmQLFTest()`, limma `eBayes()`

9. **Pathway/Gene Set Enrichment Analysis**
   - **Biological intent:** Identify enriched biological pathways
   - **Common methods:**
     - Over-representation analysis (ORA)
     - Gene Set Enrichment Analysis (GSEA)
     - GO enrichment
   - **Tools:** clusterProfiler, GSEA, GOseq, enrichR

10. **Visualization**
    - **Common plots:**
      - MA plot (log fold-change vs mean expression)
      - Volcano plot (log fold-change vs -log10 p-value)
      - Heatmap (top DE genes)
      - PCA (sample relationships)

### Typical Parameter Ranges

| Step | Parameter | Typical Range | Notes |
|------|-----------|---------------|-------|
| Filtering | min_counts | 10-20 | Total across samples |
| Filtering | min_cpm | 1-2 | In at least N samples |
| DE | padj_threshold | 0.01-0.05 | FDR cutoff |
| DE | lfc_threshold | 1-2 | Log2 fold-change |
| GSEA | fdr_threshold | 0.05-0.25 | Pathway significance |

## Spatial Transcriptomics (Visium) Workflow

### Standard Pipeline Steps

1. **Space Ranger Execution**
   - **Biological intent:** Process raw sequencing data and align to tissue image
   - **Key steps:**
     - Fiducial detection and alignment
     - Tissue detection
     - Read mapping and counting
   - **Output:** `outs/` directory with matrices and spatial metadata

2. **Quality Control**
   - **Biological intent:** Assess data quality and tissue detection
   - **Common checks:**
     - Fiducial alignment quality (check `aligned_fiducials.jpg`)
     - Tissue detection accuracy (check `detected_tissue_image.jpg`)
     - Median genes per spot
     - Median UMI per spot
   - **Tools:** Space Ranger `web_summary.html`, Seurat/Scanpy QC functions

3. **Load Data**
   - **Biological intent:** Load expression matrix with spatial coordinates
   - **Tools:** Seurat `Load10X_Spatial()`, Scanpy `sc.read_visium()`

4. **Spot Filtering (optional)**
   - **Biological intent:** Remove low-quality spots
   - **Common thresholds:**
     - Min genes per spot: 200-500
     - Min UMI per spot: 500-1000
     - Max mitochondrial %: 20-30%
   - **Note:** Less stringent than scRNA-seq (spots contain multiple cells)

5. **Normalization**
   - **Biological intent:** Account for sequencing depth differences
   - **Common methods:**
     - SCTransform (recommended for spatial)
     - Log-normalization
   - **Tools:** Seurat `SCTransform()`, Scanpy `sc.pp.normalize_total()` + `sc.pp.log1p()`

6. **Dimensionality Reduction & Clustering**
   - **Biological intent:** Identify spatial domains/regions
   - **Same as scRNA-seq:** PCA → UMAP → clustering
   - **Note:** Clusters often correspond to anatomical regions

7. **Spatially Variable Genes**
   - **Biological intent:** Identify genes with spatial expression patterns
   - **Common methods:**
     - Seurat: markvariogram, moransi
     - SpatialDE
     - SPARK
   - **Tools:** Seurat `FindSpatiallyVariableFeatures()`, SpatialDE, SPARK

8. **Spatial Analysis**
   - **Biological intent:** Analyze spatial relationships
   - **Common analyses:**
     - Spatial autocorrelation (Moran's I)
     - Neighborhood enrichment (cell type co-localization)
     - Spatial domains
   - **Tools:** Squidpy `sq.gr.spatial_autocorr()`, `sq.gr.nhood_enrichment()`

9. **Deconvolution (optional)**
   - **Biological intent:** Infer cell type composition of each spot
   - **Common methods:**
     - SPOTlight
     - cell2location
     - RCTD
     - Tangram
   - **Requires:** scRNA-seq reference

10. **Visualization**
    - **Common plots:**
      - Spatial feature plots (gene expression on tissue)
      - Spatial dim plots (clusters on tissue)
      - Heatmaps (spatially variable genes)

### Typical Parameter Ranges

| Step | Parameter | Typical Range | Notes |
|------|-----------|---------------|-------|
| QC | min_genes | 200-500 | Less stringent than scRNA |
| QC | max_mt% | 20-30% | Higher than scRNA |
| Spatial features | n_features | 1000-3000 | Spatially variable genes |
| Clustering | resolution | 0.3-0.8 | Often lower than scRNA |

## Xenium (In Situ Spatial) Workflow

### Standard Pipeline Steps

1. **Xenium Onboard Analysis**
   - **Biological intent:** Decode transcripts and segment cells
   - **Key steps:**
     - Transcript decoding
     - Cell segmentation (nucleus + cell boundaries)
     - Transcript-to-cell assignment
   - **Output:** `transcripts.parquet`, `cell_feature_matrix.h5`, `cells.zarr.zip`

2. **Segmentation Provenance Check**
   - **Biological intent:** Determine if segmentation was modified
   - **Check for:**
     - Xenium Ranger resegmentation
     - Imported segmentation (Cellpose, StarDist, etc.)
   - **Critical:** Segmentation affects all downstream results

3. **Quality Control**
   - **Biological intent:** Assess data quality
   - **Common checks:**
     - Total transcripts detected
     - Transcripts per cell
     - Cells detected
     - Transcript quality values (QV)
   - **Tools:** Check `metrics_summary.csv`

4. **Load Data**
   - **Biological intent:** Load cell × gene matrix with spatial coordinates
   - **Tools:** Seurat `LoadXenium()`, Squidpy `sq.read.xenium()`, SpatialData `sd.read_xenium()`

5. **Cell Filtering**
   - **Biological intent:** Remove low-quality cells
   - **Common thresholds:**
     - Min transcripts per cell: 20-100 (lower than scRNA-seq)
     - Min cells per gene: 5-10
   - **Note:** Targeted panel (300-500 genes) vs whole transcriptome

6. **Normalization**
   - **Biological intent:** Account for cell size and transcript capture differences
   - **Common methods:**
     - Log-normalization
     - SCTransform
   - **Tools:** Same as scRNA-seq

7. **Dimensionality Reduction & Clustering**
   - **Biological intent:** Identify cell types
   - **Same as scRNA-seq:** PCA → UMAP → clustering
   - **Note:** Single-cell resolution (unlike Visium spots)

8. **Spatial Neighbor Graph**
   - **Biological intent:** Define spatial relationships between cells
   - **Common parameters:**
     - k-nearest neighbors: 6-15
     - Distance threshold: 50-100 µm
   - **Tools:** Squidpy `sq.gr.spatial_neighbors()`

9. **Spatial Analysis**
   - **Biological intent:** Analyze spatial patterns at single-cell resolution
   - **Common analyses:**
     - Spatial autocorrelation (Moran's I)
     - Neighborhood enrichment (cell type co-localization)
     - Ligand-receptor interactions
     - Spatial domains
   - **Tools:** Squidpy `sq.gr.spatial_autocorr()`, `sq.gr.nhood_enrichment()`, `sq.gr.ligrec()`

10. **Subcellular Analysis (optional)**
    - **Biological intent:** Analyze transcript localization within cells
    - **Analyses:**
      - Nuclear vs cytoplasmic transcripts
      - Transcript clustering within cells
    - **Requires:** Matching transcript coordinates to nucleus boundaries

11. **Visualization**
    - **Common plots:**
      - Spatial scatter plots (cells colored by type)
      - Feature plots (gene expression on tissue)
      - Neighborhood enrichment heatmaps
      - Ligand-receptor networks

### Typical Parameter Ranges

| Step | Parameter | Typical Range | Notes |
|------|-----------|---------------|-------|
| QC | min_transcripts | 20-100 | Lower than scRNA (targeted panel) |
| QC | min_cells | 5-10 | Per gene |
| Spatial neighbors | k_neighbors | 6-15 | Depends on tissue density |
| Spatial neighbors | distance | 50-100 µm | Cell-cell distance threshold |

## Common Cross-Workflow Patterns

### Random Seeds (Reproducibility)

**Critical for reproducibility:**
- Clustering algorithms (Louvain, Leiden)
- UMAP/tSNE
- Subsampling operations

**Common patterns:**
```r
# R
set.seed(42)
seurat_obj <- RunUMAP(seurat_obj, seed.use = 42)
seurat_obj <- FindClusters(seurat_obj, random.seed = 42)
```

```python
# Python
import random
import numpy as np
random.seed(42)
np.random.seed(42)
sc.tl.umap(adata, random_state=42)
sc.tl.leiden(adata, random_state=42)
```

**Report if missing:** Note absence of random seeds as reproducibility risk.

### Batch Effects

**Common batch variables:**
- Sequencing run
- Sample preparation date
- Donor/patient
- Tissue section

**Detection:**
- PCA colored by batch
- UMAP colored by batch
- Batch-specific clustering

**Correction methods:**
- Harmony (most common)
- ComBat
- Seurat integration
- scVI

### Multi-Sample Integration

**Common scenarios:**
- Multiple patients/donors
- Multiple conditions
- Multiple time points

**Approaches:**
1. **Merge then analyze** (if no batch effects)
2. **Integrate then analyze** (if batch effects present)
3. **Analyze separately then compare** (for DE across conditions)

## Workflow Detection Strategy

When analyzing a repository:

1. **Identify dataset type** from file signatures:
   - `.h5ad`, `.rds` with cell × gene matrix → scRNA-seq
   - `.bam`, `.fastq` → Bulk RNA-seq
   - Space Ranger `outs/` → Visium
   - Xenium outputs → Xenium

2. **Map steps to expected workflow** from this reference

3. **Extract parameters** for each step

4. **Note deviations** from standard workflow (may indicate novel methods or errors)

5. **Check for missing steps** (e.g., no batch correction when multiple batches present)
