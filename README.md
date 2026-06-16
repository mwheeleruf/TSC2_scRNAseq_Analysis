# TSC2_scRNAseq_Analysis

Code for analyzing scRNAseq data from the TSC2 project.

## Seurat workflow

The repository now includes a Seurat-based R workflow at
`scripts/seurat_scRNAseq_pipeline.R`
for combining 6 samples across 2 groups and running:

- sample loading and merge
- integration (`LogNormalize` or `SCT`)
- normalization and scaling
- clustering and UMAP
- cluster naming
- cluster feature marker discovery
- two-group differential gene expression analysis

### Example usage

```r
source("scripts/seurat_scRNAseq_pipeline.R")

sample_paths <- c(
  "/path/to/sample_1",
  "/path/to/sample_2",
  "/path/to/sample_3",
  "/path/to/sample_4",
  "/path/to/sample_5",
  "/path/to/sample_6"
)

sample_names <- c("sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6")
groups <- c("group_a", "group_a", "group_a", "group_b", "group_b", "group_b")

results <- run_seurat_scRNAseq_pipeline(
  sample_paths = sample_paths,
  sample_names = sample_names,
  groups = groups,
  feature_genes = c("EPCAM", "COL1A1"),
  normalization_method = "LogNormalize",
  project_name = "TSC2_scRNAseq",
  output_dir = "results/seurat_pipeline"
)
```

The workflow writes merged and integrated Seurat objects, marker tables, group
differential expression results, and UMAP/feature plot PDFs into the selected
output directory. If you want cell-type labels instead of the default
`Cluster_<id>` names, pass a named `cluster_names` vector after you inspect the
identified cluster IDs.
