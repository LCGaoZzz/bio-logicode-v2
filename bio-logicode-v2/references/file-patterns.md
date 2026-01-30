# File Pattern Reference for Bio-Logicode

## File Types to Skip (Never Read)

These are binary or large data files that should never be loaded into context:

### Notebooks and Data Formats
- `*.ipynb` - Jupyter notebooks (use text scripts instead)
- `*.h5ad` - AnnData HDF5 files
- `*.h5` - HDF5 data files
- `*.rds` - R data serialization
- `*.loom` - Loom data format

### Documents and Images
- `*.pdf` - PDF documents
- `*.png`, `*.jpg`, `*.jpeg` - Image files

### Archives
- `*.zip`, `*.tar`, `*.gz` - Compressed archives

### Genomics Data
- `*.bam`, `*.sam`, `*.cram` - Alignment files
- `*.fastq`, `*.fq` - Sequencing reads
- `*.fasta`, `*.fa` - Reference sequences
- `*.mtx` - Matrix market format

### Model Files
- `*.pt`, `*.pth` - PyTorch model weights

### Spatial Data Files
- `*.parquet` - Parquet files (e.g., Xenium transcripts, cells)
- `*.zarr`, `*.zarr.zip` - Zarr archives (e.g., Xenium cells)
- `*.ome.tif`, `*.ome.tiff` - OME-TIFF images (e.g., Xenium morphology)
- `*.cloupe` - Loupe Browser files (Space Ranger output)

## High-Priority Files (Read First)

### Tier 1: Entry Points and Documentation
1. `README.md`, `README.txt`, `README`
2. `article-brief.txt` (always read first)
3. Pipeline orchestration:
   - `Snakefile`
   - `Makefile`
   - `nextflow.config`, `*.nf`
   - `main.py`, `main.R`, `main.sh`
   - `run*.sh`, `run*.py`, `run*.R`
   - `workflow*`

### Tier 2: Configuration Files
- `*.yaml`, `*.yml`
- `*.json`
- `config.*`
- `params.*`

### Tier 3: Scripts Referenced by Entry Points
After reading entry points, follow the script call chain to find the actual implementation files.

## File Naming Patterns for Analysis Steps

Common naming conventions in bioinformatics workflows:

### Quality Control
- `*qc*`, `*quality*`, `*filter*`
- `01_*`, `step1_*` (often the first step)

### Normalization
- `*norm*`, `*normalize*`, `*scale*`
- `02_*`, `step2_*`

### Dimensionality Reduction
- `*pca*`, `*umap*`, `*tsne*`
- `*dim*`, `*reduction*`

### Clustering
- `*cluster*`, `*louvain*`, `*leiden*`
- `*community*`

### Differential Expression
- `*de*`, `*deg*`, `*differential*`
- `*deseq*`, `*edger*`, `*limma*`

### Annotation
- `*annot*`, `*label*`, `*celltype*`
- `*marker*`

### Survival Analysis
- `*survival*`, `*kaplan*`, `*cox*`

### Trajectory/Pseudotime
- `*trajectory*`, `*pseudotime*`, `*monocle*`

### Spatial Analysis
- `*spatial*`, `*visium*`, `*xenium*`
- `*neighbor*`, `*nhood*`
- `*ligand*`, `*receptor*`, `*ligrec*`
- `*moran*`, `*autocorr*`

## Spatial Transcriptomics File Signatures

### Space Ranger (Visium) Output Signatures

**If you see this directory structure, Space Ranger was executed:**

```
outs/
├── web_summary.html
├── metrics_summary.csv
├── cloupe.cloupe
├── filtered_feature_bc_matrix.h5
├── filtered_feature_bc_matrix/
├── raw_feature_bc_matrix.h5
├── raw_feature_bc_matrix/
├── spatial/
│   ├── tissue_hires_image.png
│   ├── tissue_lowres_image.png
│   ├── aligned_fiducials.jpg
│   ├── detected_tissue_image.jpg
│   ├── scalefactors_json.json
│   └── tissue_positions.csv
└── analysis/
```

**Key files to check:**
1. `web_summary.html` - QC report (always check first)
2. `cloupe.cloupe` - Confirms Space Ranger completion
3. `spatial/aligned_fiducials.jpg` - Fiducial alignment QC
4. `spatial/detected_tissue_image.jpg` - Tissue detection QC
5. `filtered_feature_bc_matrix.h5` - Expression matrix (spots under tissue)

**Inference:** If these files exist, Space Ranger pipeline was executed successfully.

### Xenium Output Signatures

**If you see these files, Xenium Onboard Analysis was executed:**

```
xenium_output/
├── transcripts.parquet (or .zarr.zip)
├── cell_feature_matrix.h5
├── cell_feature_matrix/
├── cells.zarr.zip (or .parquet)
├── metrics_summary.csv
├── morphology_focus/
│   ├── morphology_focus_0000.ome.tif
│   ├── morphology_focus_0001.ome.tif
│   └── ...
└── analysis/
```

**Key files to check:**
1. `metrics_summary.csv` - QC metrics (always check first)
2. `transcripts.parquet` - Decoded transcript coordinates
3. `cell_feature_matrix.h5` - Cell × gene expression matrix
4. `cells.zarr.zip` - Cell segmentation boundaries
5. `morphology_focus/*.ome.tif` - Morphology images (DAPI, membrane)

**Inference:** If these files exist, Xenium pipeline was executed successfully.

**Segmentation provenance check:**
- If `cells.zarr.zip` timestamp > `transcripts.parquet` timestamp → Resegmentation occurred
- Look for `xeniumranger resegment` or `xeniumranger import-segmentation` commands in scripts

### Visium HD Signatures

**Additional patterns for Visium HD:**
```
outs/
├── square_002um/
│   └── filtered_feature_bc_matrix/
├── square_008um/
│   └── filtered_feature_bc_matrix/
└── square_016um/
    └── filtered_feature_bc_matrix/
```

**Inference:** Multiple bin sizes indicate Visium HD (higher resolution than standard Visium).

## Search Strategy

When looking for a specific analysis step:

1. **Check entry point first** - See if it references the step
2. **Search by number** - Look for `01_`, `02_`, etc. if numbered
3. **Search by keyword** - Use Glob with patterns like `*cluster*`, `*norm*`, `*spatial*`
4. **Check subdirectories** - Look in `scripts/`, `src/`, `analysis/`, `code/`

## Spatial Data Detection Strategy

### When Given a Directory

1. **Check for Space Ranger outputs:**
   - Look for `outs/web_summary.html` or `outs/cloupe.cloupe`
   - If found → Visium dataset, Space Ranger was executed

2. **Check for Xenium outputs:**
   - Look for `transcripts.parquet` or `cell_feature_matrix.h5` + `cells.zarr.zip`
   - If found → Xenium dataset, Xenium Onboard Analysis was executed

3. **Check for downstream analysis scripts:**
   - Look for `Load10X_Spatial()`, `LoadXenium()` (R/Seurat)
   - Look for `sc.read_visium()`, `sq.read.xenium()` (Python)
   - These indicate spatial analysis was performed

4. **Infer pipeline stages from artifacts:**
   - Even without source code, file signatures reveal pipeline execution
   - Document what was run based on output artifacts
   - Note in report: "Pipeline stages inferred from output artifacts"

## File Priority for Spatial Analysis

### Visium Projects

**Tier 1 (Read first):**
1. `outs/web_summary.html` - QC overview
2. `outs/metrics_summary.csv` - Machine-readable QC
3. Scripts calling `spaceranger count`
4. Scripts loading spatial data (`Load10X_Spatial()`, `read_visium()`)

**Tier 2:**
5. Normalization scripts
6. Clustering scripts
7. Spatial analysis scripts (`FindSpatiallyVariableFeatures()`, `spatial_autocorr()`)

**Tier 3:**
8. Visualization scripts
9. Deconvolution scripts (if present)

### Xenium Projects

**Tier 1 (Read first):**
1. `metrics_summary.csv` - QC overview
2. Scripts calling `xeniumranger` (resegment, import-segmentation)
3. Scripts loading Xenium data (`LoadXenium()`, `sq.read.xenium()`)

**Tier 2:**
4. Segmentation scripts (Cellpose, StarDist, custom)
5. Normalization and clustering scripts
6. Spatial neighbor graph construction

**Tier 3:**
7. Spatial analysis scripts (neighborhood enrichment, ligand-receptor)
8. Visualization scripts
