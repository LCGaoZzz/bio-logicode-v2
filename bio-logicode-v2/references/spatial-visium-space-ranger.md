# Visium / Space Ranger Spatial Transcriptomics Reference

## Overview

Visium is 10x Genomics' spot-based spatial transcriptomics platform. Each spot (~55 µm diameter) captures transcripts from multiple cells (~1-10 cells per spot). Space Ranger is the analysis pipeline for Visium data.

## Space Ranger Output Signatures

### Standard Space Ranger `outs/` Directory Structure

When you encounter this directory structure, you can confidently infer Space Ranger was executed:

```
outs/
├── web_summary.html          # QC report (HTML)
├── metrics_summary.csv       # QC metrics (CSV)
├── cloupe.cloupe            # Loupe Browser file
├── filtered_feature_bc_matrix.h5      # Filtered spot × gene matrix (HDF5)
├── filtered_feature_bc_matrix/        # Filtered matrix (MEX format)
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── raw_feature_bc_matrix.h5           # Raw spot × gene matrix (HDF5)
├── raw_feature_bc_matrix/             # Raw matrix (MEX format)
├── spatial/                           # Spatial metadata
│   ├── tissue_hires_image.png        # High-res tissue image
│   ├── tissue_lowres_image.png       # Low-res tissue image
│   ├── aligned_fiducials.jpg         # Fiducial alignment QC
│   ├── detected_tissue_image.jpg     # Tissue detection QC
│   ├── scalefactors_json.json        # Image scale factors
│   └── tissue_positions.csv (or .parquet)  # Spot coordinates
├── analysis/                          # Secondary analysis
│   ├── clustering/
│   ├── diffexp/
│   ├── pca/
│   └── umap/
└── molecule_info.h5                   # Molecule-level data
```

### Key File Signatures

#### QC and Summary Files

**`web_summary.html`**
- Interactive HTML QC report
- **Always check this first** when analyzing Space Ranger outputs
- Contains: sequencing metrics, mapping metrics, tissue detection, spot counts

**`metrics_summary.csv`**
- Machine-readable QC metrics
- Key metrics: total spots, spots under tissue, median genes per spot, median UMI per spot

**`cloupe.cloupe`**
- Loupe Browser visualization file
- Contains: expression matrix, spatial coordinates, secondary analysis results
- **Presence indicates Space Ranger completed successfully**

#### Expression Matrices

**Filtered vs Raw:**
- **`filtered_feature_bc_matrix`**: Spots under tissue only (recommended for analysis)
- **`raw_feature_bc_matrix`**: All spots including background

**Matrix formats:**
- `.h5` - HDF5 format (single file, faster loading)
- MEX directory - Three files (barcodes, features, matrix) - standard 10x format

**Loading in code:**
```r
# R / Seurat
library(Seurat)
visium_data <- Load10X_Spatial("path/to/outs/")
# OR
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")

# Python / Scanpy
import scanpy as sc
adata = sc.read_visium("path/to/outs/")
# OR
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
```

#### Spatial Metadata (`spatial/` folder)

**`tissue_hires_image.png`** and **`tissue_lowres_image.png`**
- H&E stained tissue images
- High-res: ~2000×2000 pixels
- Low-res: ~600×600 pixels (used for visualization)

**`aligned_fiducials.jpg`**
- QC image showing fiducial frame alignment
- **Critical QC checkpoint:** Fiducials must be correctly detected for accurate spatial mapping
- If alignment fails, spatial coordinates will be incorrect

**`detected_tissue_image.jpg`**
- QC image showing which spots were classified as "under tissue"
- **Critical QC checkpoint:** Tissue detection determines which spots are in filtered matrix

**`scalefactors_json.json`**
- Scale factors for converting between image coordinates and spot coordinates
- Example:
  ```json
  {
    "tissue_hires_scalef": 0.17011142,
    "tissue_lowres_scalef": 0.051033426,
    "fiducial_diameter_fullres": 144.0,
    "spot_diameter_fullres": 89.0
  }
  ```

**`tissue_positions.csv`** (or `.parquet` in newer versions)
- Spot barcodes and coordinates
- Columns: barcode, in_tissue (0/1), array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
- **Key column:** `in_tissue` - determines filtered vs raw matrix

### Secondary Analysis (`analysis/` folder)

Space Ranger automatically performs secondary analysis on filtered matrix:

**`analysis/clustering/`**
- Graph-based clustering (similar to Cell Ranger)
- Files: `graphclust/clusters.csv`

**`analysis/diffexp/`**
- Differential expression between clusters
- Files: `graphclust/differential_expression.csv`

**`analysis/pca/`**
- Principal component analysis
- Files: `projection.csv`, `components.csv`, `variance.csv`

**`analysis/umap/`**
- UMAP dimensionality reduction
- Files: `projection.csv`

**Note:** These results are also embedded in `cloupe.cloupe` for Loupe Browser visualization.

## Space Ranger Pipeline Detection

### Command-Line Invocation

**Look for `spaceranger count` commands:**

```bash
# Standard invocation
spaceranger count \
  --id=sample_id \
  --transcriptome=/path/to/refdata \
  --fastqs=/path/to/fastqs \
  --sample=sample_name \
  --image=/path/to/image.tif \
  --slide=slide_serial \
  --area=capture_area \
  --localcores=8 \
  --localmem=64
```

**Key parameters to extract:**
- `--id`: Sample identifier (becomes output directory name)
- `--transcriptome`: Reference genome used
- `--image`: H&E image file
- `--slide` and `--area`: Slide serial number and capture area (for spatial alignment)

### Configuration Files

**Look for:**
- `spaceranger.yaml` or similar config files
- Pipeline orchestration files (Snakemake, Nextflow) calling Space Ranger

### Log Files

**Space Ranger produces:**
- `_log` file in output directory
- `_cmdline` file with exact command used

## Matrix Choice: Filtered vs Raw

**Critical decision point:** Which matrix was used for downstream analysis?

### Detection Strategy

**If you see:**
```r
# R
Load10X_Spatial("outs/")  # Uses filtered by default
Read10X("outs/filtered_feature_bc_matrix/")  # Explicit filtered
Read10X("outs/raw_feature_bc_matrix/")  # Explicit raw
```

```python
# Python
sc.read_visium("outs/")  # Uses filtered by default
sc.read_10x_h5("filtered_feature_bc_matrix.h5")  # Explicit filtered
sc.read_10x_h5("raw_feature_bc_matrix.h5")  # Explicit raw
```

**Report format:**
```markdown
**Matrix choice:**
- Matrix used: Filtered (spots under tissue only)
- Evidence: [Script: analysis/load_data.R | Lines: 10-12]
- Spot count: [from metrics_summary.csv or code output]
```

## QC Checkpoints in Space Ranger Workflows

### Critical QC Steps to Document

1. **Fiducial alignment**
   - Check: `spatial/aligned_fiducials.jpg`
   - Look for: Comments about alignment quality, manual alignment steps

2. **Tissue detection**
   - Check: `spatial/detected_tissue_image.jpg`
   - Look for: Manual tissue selection, tissue detection parameters

3. **Sequencing metrics**
   - Check: `web_summary.html` or `metrics_summary.csv`
   - Key metrics: reads mapped to transcriptome, median genes per spot, median UMI per spot

4. **Spot filtering**
   - Look for: Additional filtering beyond Space Ranger's tissue detection
   - Common filters: min genes per spot, min UMI per spot, mitochondrial percentage

### Example QC Code Patterns

```r
# R / Seurat
visium_obj <- Load10X_Spatial("outs/")

# QC filtering
visium_obj <- subset(visium_obj,
                     subset = nFeature_Spatial > 200 &
                              nCount_Spatial > 500 &
                              percent.mt < 20)
```

```python
# Python / Scanpy
adata = sc.read_visium("outs/")

# QC filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=500)
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
```

## Visium HD Nuances (Optional)

Visium HD provides higher resolution with smaller spots (2 µm or 8 µm bins). Space Ranger outputs differ:

### Visium HD Output Differences

**Binned outputs:**
- Multiple resolutions: `square_002um/`, `square_008um/`, `square_016um/`
- Each bin size has its own `filtered_feature_bc_matrix/`

**Segmented outputs (if cell segmentation performed):**
- `cell_feature_matrix/` - Cell-level expression matrix
- Requires additional segmentation step

**Detection:**
- Look for multiple bin directories in `outs/`
- Check for `--bin-size` parameter in `spaceranger count` command

**Report format:**
```markdown
**Visium HD detected:**
- Bin sizes available: 2µm, 8µm, 16µm
- Bin size used for analysis: 8µm
- Evidence: [Script: analysis/load_data.R | Lines: 10-12]
```

## Common Visium Analysis Workflow

**Typical pipeline:**

1. **Space Ranger execution** (`spaceranger count`)
2. **Load data** (Seurat `Load10X_Spatial()` or Scanpy `read_visium()`)
3. **QC filtering** (min genes, min UMI, mitochondrial %)
4. **Normalization** (SCTransform, log-normalization)
5. **Dimensionality reduction** (PCA, UMAP)
6. **Clustering** (Louvain, Leiden)
7. **Spatial analysis:**
   - Spatially variable genes (`FindSpatiallyVariableFeatures()`, `SpatialDE`)
   - Spatial autocorrelation (Moran's I)
   - Neighborhood analysis
8. **Visualization** (spatial feature plots, spatial dim plots)

## Visium Tool Ecosystem

### R / Seurat

**Loading:**
```r
library(Seurat)
visium_obj <- Load10X_Spatial(
  data.dir = "outs/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE
)
```

**Key Seurat spatial functions:**
- `SpatialFeaturePlot()` - Visualize gene expression on tissue
- `SpatialDimPlot()` - Visualize clusters on tissue
- `FindSpatiallyVariableFeatures()` - Identify spatially variable genes

### Python / Scanpy / Squidpy

**Loading:**
```python
import scanpy as sc
import squidpy as sq

# Scanpy (basic)
adata = sc.read_visium("outs/")

# Squidpy (enhanced spatial)
adata = sq.read.visium("outs/")
```

**Key spatial analysis functions:**
- `sq.gr.spatial_neighbors()` - Build spatial neighbor graph
- `sq.gr.spatial_autocorr()` - Moran's I, Geary's C
- `sq.gr.nhood_enrichment()` - Spot neighborhood enrichment
- `sq.pl.spatial_scatter()` - Spatial visualization

## Detection Checklist for Reports

When analyzing Visium data, always document:

- [ ] **Space Ranger execution confirmed:** Evidence of `spaceranger count` command or `outs/` directory
- [ ] **Matrix choice:** Filtered or raw? Which was used for analysis?
- [ ] **QC checkpoints:** Fiducial alignment, tissue detection quality
- [ ] **Spot filtering:** Additional filtering beyond Space Ranger's tissue detection?
- [ ] **Normalization method:** SCTransform, log-normalization, other?
- [ ] **Spatial analysis:** Spatially variable genes, spatial autocorrelation, neighborhood analysis?
- [ ] **Visualization tool:** Seurat, Scanpy, Squidpy, Loupe Browser?

## Example Report Section

```markdown
### Visium Spatial Transcriptomics Processing

**Pipeline:** Space Ranger v2.0.0

**Evidence:**
- [Artifact: outs/web_summary.html - Space Ranger execution confirmed]
- [Artifact: outs/cloupe.cloupe - Secondary analysis completed]
- [Script: scripts/run_spaceranger.sh | Lines: 5-15]

**Space Ranger parameters:**
- Reference: GRCh38-2020-A
- Slide: V19L29-097
- Capture area: A1
- Evidence: [Script: scripts/run_spaceranger.sh | Lines: 8-10]

**QC checkpoints:**
- Fiducial alignment: Passed (visual inspection of aligned_fiducials.jpg)
- Tissue detection: 2,847 spots under tissue (out of 4,992 total)
- Median genes per spot: 3,245
- Median UMI per spot: 8,932
- Evidence: [Artifact: outs/metrics_summary.csv]

**Matrix choice:**
- Matrix used: Filtered (spots under tissue only)
- Evidence: [Script: analysis/load_data.R | Lines: 10-12]
- Additional filtering: min_genes=200, min_UMI=500, max_mt=20%
- Evidence: [Script: analysis/qc_filter.R | Lines: 25-30]

**Downstream analysis:**
- Tool: Seurat v4.3.0 (R)
- Loading method: `Load10X_Spatial()`
- Normalization: SCTransform
- Spatial analysis: FindSpatiallyVariableFeatures (method=markvariogram)
- Evidence: [Script: analysis/spatial_analysis.R | Lines: 40-85]

**Key parameters:**
- Normalization: SCTransform (default parameters)
- Clustering resolution: 0.8
- Spatially variable genes: top 2000 features
- Evidence: [Script: analysis/spatial_analysis.R | Lines: 50-65]
```

## Key Differences from Xenium

| Aspect | Visium (Spot-based) | Xenium (In situ) |
|--------|---------------------|------------------|
| Resolution | ~55 µm spots (~10 cells) | Single-cell / subcellular |
| Transcripts | Whole transcriptome (unbiased) | Targeted panel (~300-500 genes) |
| Output | Spot × gene matrix | Cell × gene matrix + transcript coordinates |
| Segmentation | Not required (spots are units) | Critical (defines cells) |
| Spatial analysis | Spot-level | Single-cell level |
| Key files | `filtered_feature_bc_matrix/`, `spatial/` | `transcripts.parquet`, `cells.zarr.zip` |
| Pipeline | Space Ranger | Xenium Onboard Analysis + Xenium Ranger |
