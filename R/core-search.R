#' Derive BCMP k-search bounds from batch sizes.
#'
#' Computes search bounds from batch sizes without depending on a single-cell
#' object ecosystem.
#'
#' @param batch_labels One non-missing, non-empty batch label per active cell.
#' @param n_total Number of active cells.
#' @param k_min Requested lower bound.
#' @param k_max Optional requested upper bound.
#' @param batch_frac_threshold_for_k_max Strict global cell-fraction threshold
#'   for batches eligible to define the automatic upper bound.
#'
#' @return A named list containing S, k_low, and k_high.
#' @noRd
derive_k_search_bounds <- function(
  batch_labels,
  n_total,
  k_min = 3L,
  k_max = NULL,
  batch_frac_threshold_for_k_max = 0.01
) {
  k_min <- validate_bcmp_integer(k_min, "k_min", minimum = 1)
  k_max <- validate_bcmp_integer(k_max, "k_max", minimum = 1, allow_null = TRUE)
  if (!is.atomic(batch_labels)) {
    stop("'batch_labels' must be an atomic vector.", call. = FALSE)
  }
  batch_labels <- as.character(batch_labels)
  if (!length(batch_labels) || anyNA(batch_labels) || any(!nzchar(trimws(batch_labels)))) {
    stop("Cannot derive k bounds: no non-missing batch counts found", call. = FALSE)
  }
  if (!is.numeric(n_total) || length(n_total) != 1L || is.na(n_total) ||
      !is.finite(n_total) || n_total < 1) {
    stop("'n_total' must be one finite positive number.", call. = FALSE)
  }
  n_total <- as.integer(n_total)

  if (!is.numeric(batch_frac_threshold_for_k_max) ||
      length(batch_frac_threshold_for_k_max) != 1L ||
      is.na(batch_frac_threshold_for_k_max) ||
      !is.finite(batch_frac_threshold_for_k_max) ||
      batch_frac_threshold_for_k_max < 0 ||
      batch_frac_threshold_for_k_max >= 1) {
    stop(
      "'batch_frac_threshold_for_k_max' must be a finite numeric fraction in [0, 1)",
      call. = FALSE
    )
  }

  batch_sizes <- table(batch_labels)
  k_low <- max(as.integer(k_min), 2L)
  if (is.null(k_max)) {
    eligible_batch_sizes <- batch_sizes[
      batch_sizes > n_total * batch_frac_threshold_for_k_max
    ]
    if (!length(eligible_batch_sizes)) {
      stop(
        sprintf(
          "Cannot derive k bounds: no batches exceed %.2f%% of total cells",
          100 * batch_frac_threshold_for_k_max
        ),
        call. = FALSE
      )
    }
    min_batch_n <- max(as.integer(min(eligible_batch_sizes)), 1L)
    k_high <- min(min_batch_n - 1L, n_total - 1L)
  } else {
    k_high <- min(as.integer(k_max), n_total - 1L)
  }

  list(
    S = as.integer(length(batch_sizes)),
    k_low = as.integer(k_low),
    k_high = as.integer(k_high)
  )
}
