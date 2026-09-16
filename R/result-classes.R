#' Stable fields in a BCMP final selection.
#'
#' @noRd
BCMP_SELECTION_FIELDS <- c(
  "selected_k",
  "n_domains",
  "mixing_status",
  "n_evaluated_domains",
  "n_residual_domains",
  "n_residual_cells",
  "residual_cell_frac"
)

#' Stable columns in a BCMP search trace.
#'
#' @noRd
BCMP_SEARCH_TRACE_FIELDS <- c(
  "iter_id",
  "phase",
  "k",
  "mixing_status",
  "n_domains",
  "n_failed_domains",
  "n_residual_domains",
  "actual_residual_cell_frac"
)

#' Construct a final-selection record.
#'
#' @param selected_k Representation scale selected by the search.
#' @param n_domains Number of domains in the selected partition.
#' @param mixing_status Selected partition's mixing outcome.
#' @param n_evaluated_domains Number of domains used in the mixing decision.
#' @param n_residual_domains Number of residual domains.
#' @param n_residual_cells Number of active cells in residual domains.
#' @param residual_cell_frac Fraction of active cells in residual domains.
#'
#' @return A named list with the stable selection field order.
#' @noRd
new_bcmp_selection <- function(
  selected_k,
  n_domains,
  mixing_status,
  n_evaluated_domains,
  n_residual_domains,
  n_residual_cells,
  residual_cell_frac
) {
  list(
    selected_k = selected_k,
    n_domains = n_domains,
    mixing_status = mixing_status,
    n_evaluated_domains = n_evaluated_domains,
    n_residual_domains = n_residual_domains,
    n_residual_cells = n_residual_cells,
    residual_cell_frac = residual_cell_frac
  )
}

#' Validate a standard result's shared fields.
#'
#' @param selection A final-selection record.
#' @param search_trace A search-trace data frame.
#'
#' @return `NULL`, invisibly.
#' @noRd
validate_result_contract <- function(selection, search_trace) {
  if (!is.list(selection) || !identical(names(selection), BCMP_SELECTION_FIELDS)) {
    stop("`selection` must use the stable BCMP selection field order.", call. = FALSE)
  }
  if (!is.data.frame(search_trace) ||
      !identical(names(search_trace), BCMP_SEARCH_TRACE_FIELDS)) {
    stop("`search_trace` must use the stable BCMP search-trace column order.", call. = FALSE)
  }
  invisible(NULL)
}

#' Construct a standard embedding result.
#'
#' This internal constructor preserves the documented BCMP result field order.
#' It does not run the BCMP algorithm.
#'
#' @param labels A cell-ID-named character vector of BCMP domain labels.
#' @param selection A named list describing the selected partition.
#' @param search_trace A data frame with one row per evaluated `k`.
#'
#' @return An object of class `EmbeddingResult`.
#' @noRd
new_embedding_result <- function(labels, selection, search_trace) {
  validate_result_contract(selection, search_trace)
  structure(
    list(
      labels = labels,
      selection = selection,
      search_trace = search_trace
    ),
    class = "EmbeddingResult"
  )
}

#' Construct a debug embedding result.
#'
#' @inheritParams new_embedding_result
#' @param debug A named list of candidate-level diagnostics.
#'
#' @return An object with classes `EmbeddingDebugResult` and
#'   `EmbeddingResult`.
#' @noRd
new_embedding_debug_result <- function(labels, selection, search_trace, debug) {
  validate_result_contract(selection, search_trace)
  structure(
    list(
      labels = labels,
      selection = selection,
      search_trace = search_trace,
      debug = debug
    ),
    class = c("EmbeddingDebugResult", "EmbeddingResult")
  )
}

#' Construct a standard object-workflow result.
#'
#' @param object A copied supported single-cell object.
#' @param selection A named list describing the selected partition.
#' @param search_trace A data frame with one row per evaluated `k`.
#'
#' @return An object of class `ObjectResult`.
#' @noRd
new_object_result <- function(object, selection, search_trace) {
  validate_result_contract(selection, search_trace)
  structure(
    list(
      object = object,
      selection = selection,
      search_trace = search_trace
    ),
    class = "ObjectResult"
  )
}

#' Construct a debug object-workflow result.
#'
#' @inheritParams new_object_result
#' @param debug A named list of candidate-level diagnostics.
#'
#' @return An object with classes `ObjectDebugResult` and `ObjectResult`.
#' @noRd
new_object_debug_result <- function(object, selection, search_trace, debug) {
  validate_result_contract(selection, search_trace)
  structure(
    list(
      object = object,
      selection = selection,
      search_trace = search_trace,
      debug = debug
    ),
    class = c("ObjectDebugResult", "ObjectResult")
  )
}
