BCMP_OBJECT_MIN_ACTIVE_CELLS <- 100L

validate_bcmp_object_ids <- function(ids, name) {
  ids <- as.character(ids)
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("`", name, "` must contain unique, non-missing, non-empty IDs.", call. = FALSE)
  }
  ids
}

validate_bcmp_object_batches <- function(batch_labels, cell_ids, batch_key) {
  if (anyNA(batch_labels)) {
    stop("Batch metadata column `", batch_key, "` must contain one non-missing, non-empty label per cell.", call. = FALSE)
  }
  batch_labels <- as.character(batch_labels)
  if (length(batch_labels) != length(cell_ids) || anyNA(batch_labels) || any(!nzchar(batch_labels))) {
    stop("Batch metadata column `", batch_key, "` must contain one non-missing, non-empty label per cell.", call. = FALSE)
  }
  stats::setNames(batch_labels, cell_ids)
}

resolve_bcmp_excluded_cells <- function(exclude_cells, cell_ids) {
  if (is.null(exclude_cells)) {
    return(character())
  }
  if (is.logical(exclude_cells)) {
    if (length(exclude_cells) != length(cell_ids) || anyNA(exclude_cells)) {
      stop("Logical `exclude_cells` must have one non-missing value per cell.", call. = FALSE)
    }
    return(cell_ids[exclude_cells])
  }
  if (!is.character(exclude_cells) && !is.factor(exclude_cells)) {
    stop("`exclude_cells` must be NULL, a character vector of cell IDs, or a logical cell mask.", call. = FALSE)
  }
  requested <- as.character(exclude_cells)
  if (anyNA(requested) || any(!nzchar(requested))) {
    stop("`exclude_cells` must not contain missing or empty cell IDs.", call. = FALSE)
  }
  unknown <- setdiff(unique(requested), cell_ids)
  if (length(unknown)) {
    preview <- paste(utils::head(unknown, 5L), collapse = ", ")
    suffix <- if (length(unknown) > 5L) ", ..." else ""
    stop("`exclude_cells` contains cells not present in the object: ", preview, suffix, call. = FALSE)
  }
  cell_ids[cell_ids %in% unique(requested)]
}

embed_bcmp_active_rows <- function(all_cell_ids, active_cell_ids, active_matrix, prefix) {
  active_matrix <- as.matrix(active_matrix)
  if (nrow(active_matrix) != length(active_cell_ids)) {
    stop("BCMP active embedding rows must match active cell IDs.", call. = FALSE)
  }
  full <- matrix(
    NA_real_,
    nrow = length(all_cell_ids),
    ncol = ncol(active_matrix),
    dimnames = list(all_cell_ids, paste0(prefix, "_", seq_len(ncol(active_matrix))))
  )
  full[match(active_cell_ids, all_cell_ids), ] <- active_matrix
  full
}

run_bcmp_object_workflow <- function(
  raw_counts,
  batch_labels,
  var_names,
  cell_ids,
  batch_key,
  exclude_cells,
  partition_n_pcs,
  partition_n_hvg,
  k_min,
  k_max,
  batch_frac_threshold_for_k_max,
  min_batch_coverage,
  min_cells_in_domain_per_batch,
  max_underrepresentation_fold,
  max_residual_cell_frac,
  seed,
  output_level,
  verbose
) {
  validate_bcmp_partition_parameters(
    k_min, k_max, min_batch_coverage, min_cells_in_domain_per_batch,
    max_underrepresentation_fold, seed
  )
  validate_bcmp_integer(partition_n_pcs, "partition_n_pcs", minimum = 1)
  validate_bcmp_integer(partition_n_hvg, "partition_n_hvg", minimum = 1)
  cell_ids <- validate_bcmp_object_ids(cell_ids, "cell_ids")
  batch_labels <- validate_bcmp_object_batches(batch_labels, cell_ids, batch_key)
  excluded_cell_ids <- resolve_bcmp_excluded_cells(exclude_cells, cell_ids)
  active_cell_ids <- cell_ids[!(cell_ids %in% excluded_cell_ids)]
  if (length(active_cell_ids) < BCMP_OBJECT_MIN_ACTIVE_CELLS) {
    stop(
      "Need at least ", BCMP_OBJECT_MIN_ACTIVE_CELLS,
      " active cells after applying `exclude_cells`; got ", length(active_cell_ids),
      call. = FALSE
    )
  }
  active_idx <- match(active_cell_ids, cell_ids)
  composed <- run_bcmp_from_counts(
    raw_counts = raw_counts[active_idx, , drop = FALSE],
    batch_labels = unname(batch_labels[active_cell_ids]),
    var_names = var_names,
    cell_ids = active_cell_ids,
    partition_n_pcs = partition_n_pcs,
    partition_n_hvg = partition_n_hvg,
    k_min = k_min,
    k_max = k_max,
    batch_frac_threshold_for_k_max = batch_frac_threshold_for_k_max,
    min_batch_coverage = min_batch_coverage,
    min_cells_in_domain_per_batch = min_cells_in_domain_per_batch,
    max_underrepresentation_fold = max_underrepresentation_fold,
    max_residual_cell_frac = max_residual_cell_frac,
    seed = seed,
    output_level = output_level,
    verbose = verbose
  )
  labels <- stats::setNames(rep(NA_character_, length(cell_ids)), cell_ids)
  labels[names(composed$core$labels)] <- as.character(composed$core$labels)
  if (length(excluded_cell_ids)) {
    labels[excluded_cell_ids] <- paste0("Excluded_", unname(batch_labels[excluded_cell_ids]))
  }
  if (anyNA(labels)) {
    stop("BCMP result labels do not cover all cells.", call. = FALSE)
  }
  workflow_debug <- NULL
  if (identical(output_level, "debug")) {
    workflow_debug <- list(
      hvg_genes = as.character(composed$preprocess$hvg_genes),
      excluded_cell_ids = as.character(excluded_cell_ids)
    )
  }
  list(
    core = composed$core,
    preprocess = composed$preprocess,
    labels = labels,
    active_cell_ids = active_cell_ids,
    excluded_cell_ids = excluded_cell_ids,
    workflow_debug = workflow_debug
  )
}

new_bcmp_object_debug <- function(core_debug, workflow_debug, parameters) {
  debug <- core_debug
  debug$parameters <- parameters
  debug$workflow <- workflow_debug
  debug
}
