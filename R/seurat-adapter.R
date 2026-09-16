clear_bcmp_seurat_reductions <- function(object) {
  reductions <- SeuratObject::Reductions(object)
  for (reduction in intersect(c("BCMP_PCA", "BCMP_UMAP"), reductions)) {
    object[[reduction]] <- NULL
  }
  object
}

write_bcmp_seurat_debug_reduction <- function(object, name, embeddings, assay, key) {
  object[[name]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, assay = assay, key = key
  )
  object
}

bcmp_seurat <- function(
  object, batch_key = "batch", assay = NULL, layer = NULL,
  partition_n_pcs = 30L, partition_n_hvg = 2000L,
  k_min = 3L, k_max = NULL, batch_frac_threshold_for_k_max = 0.01,
  min_batch_coverage = NULL, min_cells_in_domain_per_batch = 5L,
  max_underrepresentation_fold = 10L, max_residual_cell_frac = 0.05,
  exclude_cells = NULL, seed = 236L, output_level = "standard", verbose = TRUE
) {
  validate_output_level(output_level)
  assay_name <- if (is.null(assay)) SeuratObject::DefaultAssay(object) else as.character(assay)
  if (length(assay_name) != 1L || is.na(assay_name) || !nzchar(assay_name) ||
      !(assay_name %in% SeuratObject::Assays(object))) {
    stop("`assay` must name an assay present in the Seurat object.", call. = FALSE)
  }
  if (!inherits(object[[assay_name]], "Assay5")) {
    stop("BCMP supports Seurat v5 Assay5 objects only.", call. = FALSE)
  }
  layer_name <- if (is.null(layer)) "counts" else as.character(layer)
  if (length(layer_name) != 1L || is.na(layer_name) || !nzchar(layer_name) ||
      !(layer_name %in% SeuratObject::Layers(object, assay = assay_name))) {
    stop("`layer` must name an available raw-count layer in the selected Seurat assay.", call. = FALSE)
  }
  metadata <- object[[]]
  if (!(batch_key %in% names(metadata))) {
    stop("Missing batch metadata column: ", batch_key, call. = FALSE)
  }
  cell_ids <- colnames(object)
  raw_counts <- Matrix::t(SeuratObject::LayerData(object, assay = assay_name, layer = layer_name))
  workflow <- run_bcmp_object_workflow(
    raw_counts = raw_counts, batch_labels = metadata[cell_ids, batch_key],
    var_names = colnames(raw_counts), cell_ids = cell_ids, batch_key = batch_key,
    exclude_cells = exclude_cells, partition_n_pcs = partition_n_pcs,
    partition_n_hvg = partition_n_hvg, k_min = k_min, k_max = k_max,
    batch_frac_threshold_for_k_max = batch_frac_threshold_for_k_max,
    min_batch_coverage = min_batch_coverage,
    min_cells_in_domain_per_batch = min_cells_in_domain_per_batch,
    max_underrepresentation_fold = max_underrepresentation_fold,
    max_residual_cell_frac = max_residual_cell_frac, seed = seed,
    output_level = output_level, verbose = verbose
  )
  output <- clear_bcmp_seurat_reductions(object)
  output[["bcmp_domain"]] <- factor(unname(workflow$labels), levels = unique(unname(workflow$labels)))
  if (identical(output_level, "debug")) {
    output <- write_bcmp_seurat_debug_reduction(
      output, "BCMP_PCA",
      embed_bcmp_active_rows(cell_ids, workflow$active_cell_ids, workflow$preprocess$pca_embeddings, "BCMP_PCA"),
      assay_name, "BCMPPCA_"
    )
    if (!is.null(workflow$preprocess$umap_embeddings)) {
      output <- write_bcmp_seurat_debug_reduction(
        output, "BCMP_UMAP",
        embed_bcmp_active_rows(cell_ids, workflow$active_cell_ids, workflow$preprocess$umap_embeddings, "BCMP_UMAP"),
        assay_name, "BCMPUMAP_"
      )
    }
    parameters <- list(
      batch_key = batch_key, assay = assay_name, layer = layer_name,
      partition_n_pcs = as.integer(partition_n_pcs), partition_n_hvg = as.integer(partition_n_hvg),
      k_min = as.integer(k_min), k_max = if (is.null(k_max)) NULL else as.integer(k_max),
      batch_frac_threshold_for_k_max = as.numeric(batch_frac_threshold_for_k_max),
      min_batch_coverage = if (is.null(min_batch_coverage)) NULL else as.integer(min_batch_coverage),
      min_cells_in_domain_per_batch = as.integer(min_cells_in_domain_per_batch),
      max_underrepresentation_fold = as.numeric(max_underrepresentation_fold),
      max_residual_cell_frac = as.numeric(max_residual_cell_frac),
      exclude_cells = as.character(workflow$excluded_cell_ids), seed = as.integer(seed),
      output_level = output_level, verbose = isTRUE(verbose)
    )
    debug <- new_bcmp_object_debug(workflow$core$debug, workflow$workflow_debug, parameters)
    return(new_object_debug_result(output, workflow$core$selection, workflow$core$search_trace, debug))
  }
  new_object_result(output, workflow$core$selection, workflow$core$search_trace)
}
