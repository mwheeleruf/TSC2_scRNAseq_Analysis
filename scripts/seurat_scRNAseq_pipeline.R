assert_seurat_available <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required to run this workflow.", call. = FALSE)
  }
}

validate_pipeline_inputs <- function(sample_paths, sample_names, groups) {
  if (length(sample_paths) != 6) {
    stop("Exactly 6 sample paths are required.", call. = FALSE)
  }

  if (length(sample_names) != 6) {
    stop("Exactly 6 sample names are required.", call. = FALSE)
  }

  if (length(groups) != 6) {
    stop("Exactly 6 group labels are required.", call. = FALSE)
  }

  if (length(unique(groups)) != 2) {
    stop("Exactly 2 biological groups are required.", call. = FALSE)
  }

  if (anyDuplicated(sample_names)) {
    stop("Sample names must be unique.", call. = FALSE)
  }

  missing_paths <- sample_paths[!file.exists(sample_paths)]
  if (length(missing_paths) > 0) {
    stop(
      sprintf("Sample paths do not exist: %s", paste(missing_paths, collapse = ", ")),
      call. = FALSE
    )
  }
}

read_counts_matrix <- function(sample_path) {
  counts <- Seurat::Read10X(data.dir = sample_path)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  counts
}

create_sample_object <- function(sample_path,
                                 sample_name,
                                 group,
                                 min_cells = 3,
                                 min_features = 200) {
  counts <- read_counts_matrix(sample_path)

  sample_object <- Seurat::CreateSeuratObject(
    counts = counts,
    project = group,
    min.cells = min_cells,
    min.features = min_features
  )

  sample_object$sample_id <- sample_name
  sample_object$group <- group
  sample_object
}

merge_sample_objects <- function(sample_objects, sample_names, project_name = "TSC2_scRNAseq") {
  Seurat::merge(
    x = sample_objects[[1]],
    y = sample_objects[-1],
    add.cell.ids = sample_names,
    project = project_name
  )
}

integrate_sample_objects <- function(sample_objects,
                                     normalization_method = c("LogNormalize", "SCT"),
                                     dims = 1:30) {
  normalization_method <- match.arg(normalization_method)

  if (normalization_method == "SCT") {
    sample_objects <- lapply(
      sample_objects,
      function(object) Seurat::SCTransform(object, verbose = FALSE)
    )

    features <- Seurat::SelectIntegrationFeatures(
      object.list = sample_objects,
      nfeatures = 3000
    )

    sample_objects <- Seurat::PrepSCTIntegration(
      object.list = sample_objects,
      anchor.features = features,
      verbose = FALSE
    )

    anchors <- Seurat::FindIntegrationAnchors(
      object.list = sample_objects,
      normalization.method = "SCT",
      anchor.features = features,
      dims = dims,
      verbose = FALSE
    )

    integrated_object <- Seurat::IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      dims = dims,
      verbose = FALSE
    )
  } else {
    sample_objects <- lapply(sample_objects, function(object) {
      object <- Seurat::NormalizeData(object, verbose = FALSE)
      Seurat::FindVariableFeatures(object, verbose = FALSE)
    })

    features <- Seurat::SelectIntegrationFeatures(
      object.list = sample_objects,
      nfeatures = 3000
    )

    anchors <- Seurat::FindIntegrationAnchors(
      object.list = sample_objects,
      anchor.features = features,
      dims = dims,
      verbose = FALSE
    )

    integrated_object <- Seurat::IntegrateData(
      anchorset = anchors,
      dims = dims,
      verbose = FALSE
    )
  }

  integrated_object
}

run_clustering_umap <- function(integrated_object,
                                normalization_method = c("LogNormalize", "SCT"),
                                dims = 1:30,
                                resolution = 0.5) {
  normalization_method <- match.arg(normalization_method)

  Seurat::DefaultAssay(integrated_object) <- "integrated"

  if (normalization_method == "LogNormalize") {
    integrated_object <- Seurat::ScaleData(integrated_object, verbose = FALSE)
  }

  integrated_object <- Seurat::RunPCA(integrated_object, npcs = max(dims), verbose = FALSE)
  integrated_object <- Seurat::FindNeighbors(integrated_object, dims = dims, verbose = FALSE)
  integrated_object <- Seurat::FindClusters(integrated_object, resolution = resolution, verbose = FALSE)
  integrated_object <- Seurat::RunUMAP(integrated_object, dims = dims, verbose = FALSE)

  integrated_object
}

assign_cluster_names <- function(integrated_object, cluster_names = NULL) {
  cluster_ids <- levels(Seurat::Idents(integrated_object))

  if (is.null(cluster_names)) {
    cluster_names <- stats::setNames(
      paste0("Cluster_", cluster_ids),
      cluster_ids
    )
  }

  if (is.null(names(cluster_names)) || any(names(cluster_names) == "")) {
    stop("cluster_names must be a named character vector keyed by cluster id.", call. = FALSE)
  }

  missing_clusters <- setdiff(cluster_ids, names(cluster_names))
  if (length(missing_clusters) > 0) {
    stop(
      sprintf(
        "cluster_names is missing labels for clusters: %s",
        paste(missing_clusters, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  rename_arguments <- c(list(object = integrated_object), as.list(cluster_names))
  integrated_object <- do.call(Seurat::RenameIdents, rename_arguments)
  integrated_object$cluster_name <- as.character(Seurat::Idents(integrated_object))

  integrated_object
}

prepare_rna_assay <- function(integrated_object) {
  Seurat::DefaultAssay(integrated_object) <- "RNA"
  integrated_object <- Seurat::NormalizeData(integrated_object, verbose = FALSE)
  integrated_object <- Seurat::FindVariableFeatures(integrated_object, verbose = FALSE)
  integrated_object <- Seurat::ScaleData(integrated_object, verbose = FALSE)
  integrated_object
}

find_feature_markers <- function(integrated_object,
                                 only_pos = TRUE,
                                 min_pct = 0.25,
                                 logfc_threshold = 0.25) {
  Seurat::FindAllMarkers(
    object = integrated_object,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )
}

run_group_differential_expression <- function(integrated_object,
                                              groups,
                                              group_column = "group",
                                              min_pct = 0.1,
                                              logfc_threshold = 0.25) {
  group_levels <- unique(groups)

  Seurat::FindMarkers(
    object = integrated_object,
    group.by = group_column,
    ident.1 = group_levels[[1]],
    ident.2 = group_levels[[2]],
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )
}

save_umap_plots <- function(integrated_object, output_dir) {
  by_group_path <- file.path(output_dir, "umap_by_group.pdf")
  by_cluster_path <- file.path(output_dir, "umap_by_cluster.pdf")

  grDevices::pdf(by_group_path, width = 8, height = 6)
  print(Seurat::DimPlot(integrated_object, reduction = "umap", group.by = "group"))
  grDevices::dev.off()

  grDevices::pdf(by_cluster_path, width = 8, height = 6)
  print(Seurat::DimPlot(integrated_object, reduction = "umap", group.by = "cluster_name", label = TRUE))
  grDevices::dev.off()
}

save_feature_plots <- function(integrated_object, feature_genes, output_dir) {
  if (length(feature_genes) == 0) {
    return(invisible(NULL))
  }

  feature_plot_path <- file.path(output_dir, "feature_plots.pdf")
  grDevices::pdf(feature_plot_path, width = 10, height = 8)

  for (feature_gene in feature_genes) {
    print(Seurat::FeaturePlot(integrated_object, features = feature_gene))
  }

  grDevices::dev.off()
}

save_pipeline_outputs <- function(merged_object,
                                  integrated_object,
                                  cluster_markers,
                                  differential_expression,
                                  output_dir,
                                  feature_genes = character()) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(merged_object, file = file.path(output_dir, "merged_seurat_object.rds"))
  saveRDS(integrated_object, file = file.path(output_dir, "integrated_seurat_object.rds"))

  utils::write.csv(cluster_markers, file = file.path(output_dir, "cluster_feature_markers.csv"), row.names = FALSE)
  utils::write.csv(differential_expression, file = file.path(output_dir, "group_differential_expression.csv"))

  save_umap_plots(integrated_object, output_dir)
  save_feature_plots(integrated_object, feature_genes, output_dir)
}

run_tsc2_scrnaseq_pipeline <- function(sample_paths,
                                       sample_names,
                                       groups,
                                       cluster_names = NULL,
                                       feature_genes = character(),
                                       normalization_method = c("LogNormalize", "SCT"),
                                       dims = 1:30,
                                       resolution = 0.5,
                                       project_name = "TSC2_scRNAseq",
                                       output_dir = "results/seurat_pipeline") {
  assert_seurat_available()
  validate_pipeline_inputs(sample_paths, sample_names, groups)

  normalization_method <- match.arg(normalization_method)

  sample_objects <- Map(
    f = create_sample_object,
    sample_path = sample_paths,
    sample_name = sample_names,
    group = groups
  )

  merged_object <- merge_sample_objects(
    sample_objects = sample_objects,
    sample_names = sample_names,
    project_name = project_name
  )

  split_objects <- Seurat::SplitObject(merged_object, split.by = "sample_id")
  integrated_object <- integrate_sample_objects(
    sample_objects = split_objects,
    normalization_method = normalization_method,
    dims = dims
  )

  integrated_object <- run_clustering_umap(
    integrated_object = integrated_object,
    normalization_method = normalization_method,
    dims = dims,
    resolution = resolution
  )

  integrated_object <- assign_cluster_names(
    integrated_object = integrated_object,
    cluster_names = cluster_names
  )

  integrated_object <- prepare_rna_assay(integrated_object)
  cluster_markers <- find_feature_markers(integrated_object)
  differential_expression <- run_group_differential_expression(integrated_object, groups = groups)

  save_pipeline_outputs(
    merged_object = merged_object,
    integrated_object = integrated_object,
    cluster_markers = cluster_markers,
    differential_expression = differential_expression,
    output_dir = output_dir,
    feature_genes = feature_genes
  )

  list(
    merged_object = merged_object,
    integrated_object = integrated_object,
    cluster_markers = cluster_markers,
    differential_expression = differential_expression
  )
}
