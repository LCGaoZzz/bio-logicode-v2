---
name: bio-logicode-v2
description: Workflow-to-code translator for bioinformatics & spatial omics projects. Produces auditable reports mapping biological intent → computational steps → concrete code/artifact evidence. Use when verifying repos match manuscript Methods, reverse engineering analyses from scripts/configs/outputs (Space Ranger/Xenium/scRNA/bulk RNA), or creating reproducible workflow documentation with citations. Supports script-based and output-based analysis (precomputed pipelines).
---

# Bio-LogiCode v2: Bioinformatics Workflow-to-Code Translator

## Overview

A workflow-to-code translator for bioinformatics & spatial omics projects. Produces an auditable report mapping biological intent → computational steps → concrete code/artifact evidence.

## Core Capabilities

### 1. Workflow-to-Code Mapping

Given a paper Methods section / protocol / analysis notes + a repo/output folder, produce a stepwise workflow map with:

- Biological intent (what question is being answered)
- Computational intent (what algorithm/method is used)
- Concrete implementation evidence (file + line, or exact output artifact)

### 2. Evidence-Grade Code Citation

Every claim must have a pointer to:

- A script/function (with line span), OR
- A config/CLI command, OR
- A standard pipeline output artifact (e.g., Space Ranger `outs/`, Xenium output bundle)

Citation format: `[Script: path/to/file.R | Lines: 45-50]` or `[Artifact: outs/web_summary.html]`

### 3. Intelligent Exploration (Works with Precomputed Pipelines)

If given outputs rather than scripts, infer pipeline stages from **file signatures**:

- Space Ranger: `web_summary.html`, `cloupe.cloupe`, `filtered_feature_bc_matrix/`, `spatial/` folder
- Xenium: `transcripts.parquet`, `cell_feature_matrix.h5`, `cells.zarr.zip`, morphology OME-TIFF
- Standard outputs: `.h5ad`, `.rds`, `.loom` files with metadata

See [references/file-patterns.md](references/file-patterns.md) for complete signatures.

### 4. Multi-Language + Multi-Ecosystem Fluency

Support for:

- R (Seurat, DESeq2, edgeR, limma)
- Python (Scanpy, Squidpy, SpatialData, scVI)
- Shell pipelines (Snakemake, Nextflow, Bash)
- Standard pipelines (Space Ranger, Xenium Ranger, Cell Ranger)

### 5. Progressive Reporting (Checkpoint System)

Always produce partial, structured deliverables as evidence accumulates (not at the end). See Workflow Methodology below.

## Domain Knowledge Playbooks

Bio-LogiCode has specialized knowledge for:

- **scRNA-seq**: QC → normalization → dimensionality reduction → clustering → DE → annotation
- **Bulk RNA-seq**: QC → alignment → quantification → normalization → DE → pathway analysis
- **Spatial Transcriptomics (Visium / Space Ranger)**: See [references/spatial-visium-space-ranger.md](references/spatial-visium-space-ranger.md)
- **In situ single-cell (Xenium / Xenium Ranger)**: See [references/spatial-xenium.md](references/spatial-xenium.md)

For common workflow patterns, see [references/common-workflows.md](references/common-workflows.md).

For tool-specific patterns (Seurat, Scanpy, Squidpy, etc.), see [references/tool-patterns.md](references/tool-patterns.md).

## Workflow Methodology

### Three-Pass Exploration

**Pass 1: Identify entrypoints + pipeline boundaries**

1. Read `article-brief.txt` or Methods section
2. Extract ≤12 key analytical steps
3. Identify dataset type (scRNA/bulk/Visium/Xenium/mixed)
4. **OUTPUT CHECKPOINT 1**: Step list + dataset type + exploration plan

**Pass 2: Build workflow outline + evidence table**

1. Locate pipeline entrypoints (Snakefile, main scripts, config files)
2. Map each step to candidate files/artifacts
3. Read high-priority files (≤15 files in this pass)
4. **OUTPUT CHECKPOINT 2**: Workflow outline with evidence plan

**Pass 3: Deep dive into key steps + parameter provenance**

1. For each step, extract key computational nodes
2. Cite specific lines and parameters
3. Note thresholds, cutoffs, model formulas
4. **OUTPUT CHECKPOINT 3**: Evidence table populated

### Two-Stage Output Pattern

**Stage A: Structured Outline (Fast, Controllable)**

Output a **Workflow Outline** with bullet points, each step paired with:

- Expected inputs/outputs
- Where to look for evidence (files, configs, logs, artifacts)
- Key parameters to extract

**Stage B: Narrative Report (Convert Outline → Prose)**

Convert outline into a readable report (Methods-like), but keep citations. This separation makes the skill robust: planning stays lightweight; prose is generated only when the structure is stable.

## Output Format (Report Contract)

**CRITICAL: Always save the final report to a file.**

When analysis is complete, you MUST:
1. Generate the complete report following the format below
2. Save it to a file named `bioinformatics-analysis-report.md` in the analysis directory
3. Inform the user of the saved file location

If the user provides a specific output filename or location, use that instead.

The report format is:

```markdown
# Bioinformatics Analysis Implementation Report

## Executive Summary
[2-3 sentences: dataset type, main workflow, key findings]

## 1. Study Overview
[One paragraph on the workflow from article-brief.txt or Methods]

## 2. Workflow Map (Step-by-Step)

### 2.1 <Step Name>
**Biological intent:** [What question is being answered]

**Computational intent:** [What method/algorithm is used]

**Implementation:** [2-3 sentences describing the code/function]

**Evidence:**
- [Script: path/to/file.R | Lines: 45-50]
- [Script: path/to/script.py | Lines: 120-125]
OR
- [Artifact: outs/web_summary.html - Space Ranger execution confirmed]
- [Artifact: outs/filtered_feature_bc_matrix.h5 - filtered matrix used]

**Key parameters:**
- Threshold: [value]
- Method: [specific algorithm/package]
- Seed: [if applicable]

**Notes:** [Only if needed - unusual choices, missing info, etc.]

### 2.2 <Step Name>
[Repeat for each step]

## 3. Evidence Table (Step ↔ File ↔ Lines/Artifact)

| Step | File/Artifact | Lines/Type | Key Parameters |
|------|---------------|------------|----------------|
| QC filtering | scripts/qc.R | 45-50 | min_genes=200, min_cells=3 |
| Normalization | scripts/norm.py | 120-125 | method=LogNormalize |
| ... | ... | ... | ... |

## 4. Parameters & Thresholds

- **QC cutoffs:** [list]
- **Normalization method:** [method + parameters]
- **Clustering resolution:** [value]
- **DE thresholds:** [FDR, log2FC cutoffs]
- **Random seeds:** [present/absent + values]

## 5. QC Checkpoints & Failure Modes

- [Which QC steps are present]
- [What failure modes are handled]
- [What validation checks exist]

## 6. Reproducibility Notes

- **Strengths:** [What makes this reproducible]
- **Risks:** [Missing seeds, hidden defaults, version dependencies]
- **External dependencies:** [Package versions, data sources]

## 7. Discrepancies / Reviewer Notes

- [Mismatches between brief and code]
- [Unusual parameter choices]
- [Missing steps or unclear implementations]

## 8. Open Questions / Missing Evidence

- [What couldn't be verified]
- [What files would answer remaining questions]
- [What should be clarified with authors]

## 9. Summary Statistics

- Files analyzed: X of Y total scripts
- Workflow coverage: [percentage or description]
- Evidence quality: [strong/moderate/weak for each step]
```

## Efficient Workflow Principles

These principles ensure reliable completion while maximizing coverage:

### Prioritize High-Signal Files

**DO NOT scan the whole repo.** Follow this priority order:

1. **First priority:**
   - `article-brief.txt` or Methods section
   - README files
   - Pipeline entrypoints: `Snakefile`, `Makefile`, `main.*`, `run*.sh`, `workflow*`, `nextflow.config`, `*.nf`
   - Config files: `*.yaml`, `*.yml`, `*.json`

2. **Second priority:**
   - Scripts referenced by entrypoints
   - Step-specific scripts (identified by naming patterns)

3. **Third priority:**
   - Utility/helper scripts (only if needed for parameter provenance)

### Read Budget (Generous Limit)

- Read at most **20 files OR 500 KB total text**, whichever comes first
- If a file is >400 KB: Read only **first 200 lines + last 100 lines** and mark as truncated
- Track progress and aim to use the full budget for comprehensive coverage
- When reaching 18/20 files or 450/500 KB, prepare to finalize after next 1-2 files

### Skip Binaries; Focus on Scripts/Configs

**DO NOT open these file types:**
```
*.ipynb *.h5ad *.h5 *.rds *.loom *.pdf *.png *.jpg *.jpeg
*.zip *.tar *.gz *.bam *.sam *.cram *.fastq *.fq *.fasta
*.fa *.mtx *.pt *.pth *.zarr
```

If needed, mention them as inputs/outputs only (do not ingest).

### Progressive Checkpoints (Minimum Contract)

- **Checkpoint 1:** Entry points + dataset type (bulk/scRNA/ST/Xenium/Visium) + exploration plan
- **Checkpoint 2:** Workflow outline + evidence plan (what files will be read next)
- **Checkpoint 3:** Evidence table populated for core steps
- **Final:** Narrative report + reproducibility/risk notes + open questions + **SAVE TO FILE**

Output checkpoints after every 3-5 file reads. Never wait for a "perfect final answer."

**IMPORTANT:** After generating the final report, you MUST save it to a file using the Write tool. Default filename: `bioinformatics-analysis-report.md` in the analysis directory.

### Error Recovery & Resilience

**If a file read fails (permissions, encoding, corruption):**

1. Log the error: "⚠️ Failed to read [file]: [reason]"
2. Mark it in "Open Questions / Missing Evidence" section
3. **Continue with next priority file** - do not stop the analysis
4. If multiple consecutive failures (3+), switch to exploration mode

**Graceful degradation:**

- Partial file reads are acceptable (use what you got)
- Missing files should be noted but not block progress
- If a critical entrypoint is missing, infer structure from file naming patterns

### Exploration Mode (for Precomputed Pipelines)

**Trigger exploration mode when:**

- Key implementations are precomputed (loading .rds, .h5ad, .pkl files)
- No clear algorithm implementation found after 10 file reads
- Article-brief describes methods but code only shows data loading
- Given output artifacts (Space Ranger `outs/`, Xenium bundle) without source code

**Exploration strategies:**

1. **Infer from file signatures:**
   - Space Ranger outputs → Space Ranger pipeline was run
   - Xenium outputs → Xenium Onboard Analysis + possible Xenium Ranger reanalysis
   - See [references/file-patterns.md](references/file-patterns.md)

2. **Search for function definitions:**
   - Use Grep to find custom functions: `function_name\s*<-\s*function` (R)
   - Search for Python functions: `def function_name`
   - Look in utility/helper script directories

3. **Check for external dependencies:**
   - Look for library imports that might contain the implementation
   - Note package versions in requirements.txt, renv.lock, environment.yml
   - Document which external packages likely contain the algorithms

4. **Infer from comments and variable names:**
   - Read code comments for algorithm descriptions
   - Variable names like `scvi_corrected`, `kmeans_k15` reveal methods
   - Look for citations in comments (e.g., "# Using method from Smith et al. 2023")

5. **Search for configuration files:**
   - YAML/JSON configs may contain algorithm parameters
   - Look for parameter files (params.R, config.py)

**Output format for exploration findings:**

```markdown
**Implementation (exploration mode):**
The actual algorithm is precomputed in [file.rds]. Based on code context:
- Variable name `scvi_latent` suggests scVI batch correction was used
- Import statement `library(scvi)` confirms the package
- Likely parameters: [inferred from downstream usage]
- **Recommendation:** Check upstream preprocessing scripts or contact authors for implementation details
```

## Spatial Transcriptomics Support

Bio-LogiCode has specialized support for spatial omics workflows:

### Visium / Space Ranger

For spot/bin-based spatial transcriptomics. See [references/spatial-visium-space-ranger.md](references/spatial-visium-space-ranger.md) for:

- Space Ranger output signatures (`outs/` directory structure)
- QC images and fiducial alignment verification
- Matrix choice (filtered vs raw)
- Secondary analysis provenance
- Visium HD nuances

### Xenium / Xenium Ranger

For single-cell / subcellular in situ spatial transcriptomics. See [references/spatial-xenium.md](references/spatial-xenium.md) for:

- Xenium artifact signatures (transcripts, cell matrices, morphology)
- Segmentation provenance (critical for Xenium)
- Xenium Ranger pipelines (relabel, resegment, import-segmentation)
- Tool ecosystem patterns (Seurat, Squidpy, SpatialData)

### Spatial Analysis Patterns

Common downstream spatial analyses to recognize:

- **Neighbor graph construction** (spatial proximity)
- **Neighborhood enrichment** (cell type co-localization)
- **Spatial autocorrelation** (Moran's I, Geary's C)
- **Ligand-receptor interaction analysis** (cell-cell communication)

See [references/tool-patterns.md](references/tool-patterns.md) for implementation patterns.

## Key Principles

- **Always produce something** - Partial reports are better than no output
- **Checkpoint frequently** - After every 3-5 file reads
- **Use the full budget** - Aim for 18-20 files or 450-500 KB for comprehensive coverage
- **Focus on key nodes** - Not every line of code needs citation
- **Be concise** - Each step section should be 2-8 sentences max
- **Skip binaries** - Never try to read data files or compiled artifacts
- **Prioritize smartly** - Start with entrypoints and configs, not random files
- **Recover from errors** - File read failures should not stop the analysis
- **Explore when stuck** - Use exploration mode when implementations are precomputed
- **Cite precisely** - Every claim needs a file + line or artifact reference
- **SAVE THE REPORT** - Always use the Write tool to save the final report to a markdown file

## Example Workflow

**User provides:** Directory with `article-brief.txt` and `code/`

**Your process:**

1. Read `article-brief.txt` → Extract 8 steps → Identify as scRNA-seq → **OUTPUT CHECKPOINT 1**
2. Read `code/README.md` → Find entrypoint is `run_analysis.sh` → **OUTPUT CHECKPOINT 2**
3. Read `run_analysis.sh` → Identifies 3 main scripts → Build workflow outline → **OUTPUT CHECKPOINT 3**
4. Read first 3 scripts → Map steps 1-3 with citations → **OUTPUT CHECKPOINT 4**
5. Read next 3 scripts → Map steps 4-6 with citations → **OUTPUT CHECKPOINT 5**
6. Read next 3 scripts → Map steps 7-8 with citations → **OUTPUT CHECKPOINT 6**
7. **Check coverage:** If <80% of steps mapped and budget remains, continue reading
8. Generate final comprehensive report with all sections
9. **SAVE REPORT:** Use Write tool to save report as `bioinformatics-analysis-report.md` in the analysis directory
10. Inform user of the saved file location

**Spatial example (output-based):**

1. User provides Xenium output bundle (no scripts)
2. Identify `transcripts.parquet`, `cell_feature_matrix.h5`, `cells.zarr.zip` → Xenium detected → **OUTPUT CHECKPOINT 1**
3. Check for `metrics_summary.csv` → Extract QC metrics → **OUTPUT CHECKPOINT 2**
4. Infer segmentation provenance from file timestamps and structure → **OUTPUT CHECKPOINT 3**
5. Look for downstream analysis scripts (if any) → Map to Xenium outputs → **OUTPUT CHECKPOINT 4**
6. Generate report documenting pipeline stages inferred from artifacts
7. **SAVE REPORT:** Use Write tool to save report as `bioinformatics-analysis-report.md`
8. Inform user of the saved file location

## Reference Files

- [references/file-patterns.md](references/file-patterns.md) - File type patterns and prioritization
- [references/common-workflows.md](references/common-workflows.md) - scRNA, bulk RNA, ST workflow templates
- [references/tool-patterns.md](references/tool-patterns.md) - Seurat, Scanpy, Squidpy, SpatialData patterns
- [references/spatial-visium-space-ranger.md](references/spatial-visium-space-ranger.md) - Visium/Space Ranger knowledge
- [references/spatial-xenium.md](references/spatial-xenium.md) - Xenium/Xenium Ranger knowledge
