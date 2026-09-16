#' Run BCMP from an unintegrated embedding
#'
#' Run BCMP on an existing unintegrated cell-by-dimension embedding and its
#' batch labels.
#'
#' @param embedding (`matrix` / `sparseMatrix`) Numeric matrix or a numeric
#'   `Matrix::sparseMatrix` with stored numeric values (e.g. `dgCMatrix`,
#'   `dgRMatrix`, or `dgTMatrix`). Cells are rows and dimensions are columns.
#' @param batch_labels (`vector`) One non-missing, non-empty batch label per
#'   embedding row.
#' @param cell_ids (`vector` / `NULL`) Unique, non-missing cell IDs. With `NULL`, IDs
#'   are `"0"` through `"n-1"`.
#' @param k_min (`integer`) Lower bound of the representation-scale search.
#' @param k_max (`integer` / `NULL`) Optional upper bound. `NULL` derives it from
#'   batch sizes.
#' @param batch_frac_threshold_for_k_max (`numeric`) Batches at or below this global
#'   cell fraction are ignored during automatic upper-bound derivation.
#' @param min_batch_coverage (`integer` / `NULL`) Minimum number of batches required
#'   in each evaluated domain. `NULL` resolves to `max(2, ceiling(n_batches * 0.5))`.
#' @param min_cells_in_domain_per_batch (`integer`) Minimum effective cells required
#'   for a batch to count as represented in a domain.
#' @param max_underrepresentation_fold (`numeric`) Maximum allowed fold
#'   under-representation relative to global batch proportion.
#' @param max_residual_cell_frac (`numeric`) Maximum fraction of active cells that
#'   may fall in residual domains excluded from the pass/fail decision.
#' @param seed (`integer`) Random seed used by partitioning.
#' @param output_level (`character`) `"standard"` (default) or `"debug"`.
#' @param verbose (`logical`) Emit progress records when `TRUE`.
#'
#' @return `"standard"` returns an `EmbeddingResult` with `labels`, `selection`,
#'   and `search_trace`. `"debug"` returns an `EmbeddingDebugResult`, which adds
#'   `debug`.
#'
#' @export
bcmp_embedding <- function(
  embedding,
  batch_labels,
  cell_ids = NULL,
  k_min = 3L,
  k_max = NULL,
  batch_frac_threshold_for_k_max = 0.01,
  min_batch_coverage = NULL,
  min_cells_in_domain_per_batch = 5L,
  max_underrepresentation_fold = 10L,
  max_residual_cell_frac = 0.05,
  seed = 236L,
  output_level = "standard",
  verbose = TRUE
) {
  input <- validate_embedding_input(embedding, batch_labels, cell_ids)
  validate_output_level(output_level)
  run_bcmp_partition_search(
    embedding = input$embedding,
    batch_labels = input$batch_labels,
    cell_ids = input$cell_ids,
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
}

#' Validate a requested output level.
#'
#' @param output_level Requested output level.
#'
#' @return The validated output level, invisibly.
#' @noRd
validate_output_level <- function(output_level) {
  if (!is.character(output_level) || length(output_level) != 1L ||
      is.na(output_level) || !(output_level %in% c("standard", "debug"))) {
    stop("`output_level` must be either \"standard\" or \"debug\".", call. = FALSE)
  }
  invisible(output_level)
}
