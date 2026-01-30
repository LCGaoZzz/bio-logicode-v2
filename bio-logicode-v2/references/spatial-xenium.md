# Xenium In Situ Spatial Transcriptomics Reference

## Overview

Xenium is 10x Genomics' in situ single-cell/subcellular spatial transcriptomics platform. Unlike spot-based methods (Visium), Xenium provides single-molecule resolution with subcellular localization.

## Xenium Output Artifact Signatures

### Core Output Files (Xenium Onboard Analysis)

When you encounter these files, you can confidently infer a Xenium pipeline was executed:

#### Transcript Data
- **`transcripts.parquet`** or **`transcripts.zarr.zip`**
  - Decoded transcript coordinates (X, Y, Z)
  - Each row = one transcript molecule
  - Columns: transcript_id, feature_name, x_location, y_location, z_location, qv (quality value), cell_id (if assigned)
  - **Key insight:** If `cell_id` is present, transcripts have been assigned to cells via segmentation

#### Cell Feature Matrix
- **`cell_feature_matrix.h5`** - HDF5 format cell × gene matrix
- **`cell_feature_matrix/`** - MEX format directory (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
  - Standard 10x matrix format
  - Can be loaded with Seurat `Read10X()` or Scanpy `read_10x_h5()`

#### Segmentation Data
- **`cells.zarr.zip`** or **`cells.parquet`**
  - Cell and nucleus boundaries (polygons or masks)
  - Cell metadata: cell_id, x_centroid, y_centroid, nucleus_area, cell_area
  - **Critical:** Xenium produces a **flattened 2D segmentation mask** even though nucleus segmentation uses 3D information
  - Transcripts are assigned based on X/Y coordinates within cell boundaries

#### Morphology Images
- **`morphology_focus/`** directory with multi-file OME-TIFF images
  - Example files: `morphology_focus_0000.ome.tif`, `morphology_focus_0001.ome.tif`, etc.
  - DAPI nuclear staining: often in channel 0 (e.g., `ch0000_dapi.ome.tif` naming pattern)
  - Membrane staining: additional channels
  - **Format:** OME-TIFF (Open Microscopy Environment TIFF) - multi-page, multi-channel

#### QC and Metrics
- **`metrics_summary.csv`**
  - Pipeline-level QC metrics
  - Total transcripts, transcripts per cell, cells detected, etc.
  - **Always read this first** when analyzing Xenium outputs

- **`analysis/`** directory
  - Secondary analysis outputs (clustering, UMAP, etc.)
  - Similar structure to Cell Ranger secondary analysis

### File Signature Detection Logic

**If you see:**
```
transcripts.parquet (or .zarr.zip)
+ cell_feature_matrix.h5
+ cells.zarr.zip (or .parquet)
+ morphology_focus/*.ome.tif
```

**Then infer:**
- Xenium Onboard Analysis was executed
- Segmentation has been performed (cells are defined)
- Transcripts have been assigned to cells
- Morphology images are available for visualization

## Segmentation Provenance (Critical for Xenium)

Segmentation is not a side detail in Xenium—many downstream quantities depend on it. **Always include a "Segmentation Provenance" subsection in reports.**

### Three Segmentation Sources

1. **Onboard segmentation** (default)
   - Performed during Xenium Onboard Analysis
   - Uses DAPI + membrane staining
   - Produces initial `cells.zarr.zip`

2. **Xenium Ranger resegmentation** (`xeniumranger resegment`)
   - Re-runs segmentation with different parameters
   - Uses same morphology images
   - Produces new `cells.zarr.zip` and updated `cell_feature_matrix`

3. **Imported segmentation** (`xeniumranger import-segmentation`)
   - Brings in 2D nuclei/cell segmentations from community tools (Cellpose, StarDist, etc.)
   - **Recalculates all outputs that depend on segmentation:**
     - `cell_feature_matrix` (transcript-to-cell assignments)
     - `cells.zarr.zip` (cell boundaries)
     - Secondary analysis results
   - **Input requirements:** 2D segmentation masks (nucleus and/or cell)

### How to Determine Segmentation Provenance

**Check for:**

1. **Xenium Ranger logs or scripts**
   - Look for `xeniumranger resegment` or `xeniumranger import-segmentation` commands
   - Check for config files with segmentation parameters

2. **File timestamps**
   - If `cells.zarr.zip` is newer than `transcripts.parquet`, resegmentation likely occurred
   - If `cell_feature_matrix.h5` is newer than `transcripts.parquet`, reanalysis occurred

3. **Directory structure**
   - Multiple output directories (e.g., `xenium_output/`, `resegmented_output/`) suggest reanalysis

4. **Comments in scripts**
   - Look for mentions of "resegmentation", "Cellpose", "custom segmentation", etc.

**Report format:**
```markdown
**Segmentation provenance:**
- Source: [Onboard / Xenium Ranger resegment / Imported from external tool]
- Tool: [If imported: Cellpose, StarDist, custom, etc.]
- Evidence: [File timestamps, script commands, config files]
- Parameters: [If available: segmentation thresholds, model used]
```

## Xenium Ranger Pipelines

Xenium Ranger provides post-run reanalysis pipelines. **If you see `xeniumranger` commands in scripts, identify which pipeline:**

### `xeniumranger resegment`
- Re-runs segmentation with different parameters
- Uses original morphology images
- **Outputs:** New `cells.zarr.zip`, updated `cell_feature_matrix`, updated secondary analysis

### `xeniumranger import-segmentation`
- Imports 2D segmentation masks from external tools
- **Input:** Nucleus mask (required), cell mask (optional)
- **Outputs:** Same as resegment (new cells, updated matrix)
- **Use case:** When users want to use community segmentation tools (Cellpose, StarDist, etc.)

### `xeniumranger relabel`
- Relabels cell types or clusters
- Does not change segmentation or transcript assignments

### `xeniumranger rename`
- Renames features (genes) in the output
- Does not reprocess data

### Detection Strategy

**Search for:**
```bash
# In shell scripts
xeniumranger resegment --id=...
xeniumranger import-segmentation --id=... --nuclei-mask=...

# In Python scripts
import subprocess
subprocess.run(["xeniumranger", "resegment", ...])

# In config files
pipeline: xeniumranger
command: resegment
```

## Xenium Tool Ecosystem (Downstream Analysis)

After Xenium Onboard Analysis, users typically perform downstream analysis with community tools. **Recognize these patterns:**

### R / Seurat

**Loading Xenium data:**
```r
library(Seurat)

# Method 1: LoadXenium (Seurat v5+)
xenium_obj <- LoadXenium("path/to/xenium_output/", fov = "fov")

# Method 2: Load10X_Spatial (older)
xenium_obj <- Load10X_Spatial("path/to/xenium_output/")

# Method 3: Read matrix + add spatial coords
counts <- Read10X("path/to/cell_feature_matrix/")
xenium_obj <- CreateSeuratObject(counts)
# ... add spatial coordinates from cells.parquet
```

**Key Seurat spatial functions:**
- `ImageDimPlot()` - Visualize clusters on tissue
- `ImageFeaturePlot()` - Visualize gene expression on tissue
- `FindSpatiallyVariableFeatures()` - Identify spatially variable genes

### Python / Scanpy / Squidpy / SpatialData

**Loading Xenium data:**
```python
import scanpy as sc
import squidpy as sq
import spatialdata as sd

# Method 1: Scanpy (load matrix only)
adata = sc.read_10x_h5("cell_feature_matrix.h5")
# ... add spatial coordinates from cells.parquet

# Method 2: Squidpy (load matrix + spatial)
adata = sq.read.xenium("path/to/xenium_output/")

# Method 3: SpatialData (comprehensive spatial data structure)
sdata = sd.read_xenium("path/to/xenium_output/")
```

**Key spatial analysis patterns:**
- `sq.gr.spatial_neighbors()` - Build spatial neighbor graph
- `sq.gr.spatial_autocorr()` - Moran's I, Geary's C
- `sq.gr.nhood_enrichment()` - Cell type co-localization
- `sq.gr.ligrec()` - Ligand-receptor interaction analysis

### Detection Strategy

**If you see these imports/functions, expect spatial analysis:**

```r
# R patterns
library(Seurat)
LoadXenium()
ImageDimPlot()
FindSpatiallyVariableFeatures()
```

```python
# Python patterns
import squidpy as sq
import spatialdata as sd
sq.read.xenium()
sq.gr.spatial_neighbors()
sq.gr.nhood_enrichment()
sq.gr.ligrec()
```

**Report format:**
```markdown
**Downstream analysis tools:**
- Platform: [R/Python]
- Primary package: [Seurat/Scanpy/Squidpy/SpatialData]
- Loading method: [LoadXenium/Read10X/sq.read.xenium/etc.]
- Spatial analyses: [neighbor graphs, spatial autocorrelation, ligand-receptor, etc.]
- Evidence: [Script: path/to/file.R | Lines: 10-25]
```

## Xenium-Specific Analysis Patterns

### Transcript-Level Analysis

Unlike spot-based methods, Xenium provides single-molecule data. Look for:

```python
# Reading transcript data directly
import pandas as pd
transcripts = pd.read_parquet("transcripts.parquet")

# Filtering by quality
transcripts_hq = transcripts[transcripts['qv'] > 20]

# Spatial filtering (region of interest)
transcripts_roi = transcripts[
    (transcripts['x_location'] > x_min) &
    (transcripts['x_location'] < x_max)
]
```

### Subcellular Localization

Xenium provides Z-coordinates (though segmentation is 2D). Look for:

```python
# Nuclear vs cytoplasmic transcripts
# (requires matching transcript coords to nucleus boundaries)
transcripts['subcellular_location'] = assign_subcellular_location(
    transcripts, nucleus_boundaries
)
```

### High-Resolution Spatial Analysis

Xenium's single-cell resolution enables finer spatial analysis:

```python
# Cell-cell distance matrices
from scipy.spatial.distance import cdist
cell_coords = cells[['x_centroid', 'y_centroid']].values
distances = cdist(cell_coords, cell_coords)

# Neighborhood analysis at single-cell resolution
# (vs spot-level in Visium)
```

## Common Xenium Workflow Steps

**Typical analysis pipeline:**

1. **Load data** (Seurat/Scanpy/Squidpy)
2. **QC filtering** (min transcripts per cell, min cells per gene)
3. **Normalization** (log-normalization, SCTransform)
4. **Dimensionality reduction** (PCA, UMAP)
5. **Clustering** (Louvain, Leiden)
6. **Cell type annotation** (marker genes, reference mapping)
7. **Spatial analysis:**
   - Spatial neighbor graph construction
   - Spatially variable genes
   - Neighborhood enrichment (cell type co-localization)
   - Ligand-receptor interactions
8. **Visualization** (spatial plots, feature plots)

## Key Differences from Visium

| Aspect | Visium (Spot-based) | Xenium (In situ) |
|--------|---------------------|------------------|
| Resolution | ~55 µm spots (~10 cells) | Single-cell / subcellular |
| Transcripts | Whole transcriptome (unbiased) | Targeted panel (~300-500 genes) |
| Output | Spot × gene matrix | Cell × gene matrix + transcript coordinates |
| Segmentation | Not required (spots are units) | **Critical** (defines cells) |
| Spatial analysis | Spot-level | Single-cell level |
| Key files | `filtered_feature_bc_matrix/`, `spatial/` | `transcripts.parquet`, `cells.zarr.zip` |

## Checklist for Xenium Analysis Reports

When analyzing Xenium data, always address:

- [ ] **Segmentation provenance:** Onboard / resegmented / imported?
- [ ] **Segmentation tool:** If imported, which tool? (Cellpose, StarDist, custom)
- [ ] **Transcript assignment:** How were transcripts assigned to cells?
- [ ] **QC metrics:** Total transcripts, transcripts per cell, cells detected
- [ ] **Downstream tool:** Seurat / Scanpy / Squidpy / SpatialData?
- [ ] **Spatial analyses:** Neighbor graphs, spatial autocorrelation, ligand-receptor?
- [ ] **Visualization:** How are spatial plots generated?

## Example Report Section

```markdown
### Xenium Data Processing

**Dataset type:** Xenium in situ spatial transcriptomics

**Segmentation provenance:**
- Source: Imported from Cellpose
- Evidence: [Script: scripts/run_cellpose.py | Lines: 15-30]
- Xenium Ranger import: [Script: scripts/import_segmentation.sh | Lines: 5-10]
- Parameters: Cellpose model=cyto2, diameter=15

**Data loading:**
- Tool: Squidpy (Python)
- Method: `sq.read.xenium()`
- Evidence: [Script: analysis/load_data.py | Lines: 8-12]

**Spatial analysis:**
- Neighbor graph: k=6 nearest neighbors
- Neighborhood enrichment: Cell type co-localization analysis
- Ligand-receptor: CellPhoneDB database
- Evidence: [Script: analysis/spatial_analysis.py | Lines: 45-80]

**Key parameters:**
- Min transcripts per cell: 50
- Min cells per gene: 10
- Spatial neighbor k: 6
- Clustering resolution: 0.5
```
