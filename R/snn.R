#' Determine the minimum shared-neighbor count for SNN pruning.
#'
#' @param k Number of ranked neighbors per cell.
#' @param prune_snn Seurat-compatible SNN pruning threshold.
#'
#' @return An integer shared-neighbor threshold.
#' @noRd
min_shared_neighbors_for_prune <- function(k, prune_snn = 1 / 15) {
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k < 0 || k != as.integer(k)) {
    stop("'k' must be one finite non-negative integer.", call. = FALSE)
  }
  if (!is.numeric(prune_snn) || length(prune_snn) != 1L ||
      is.na(prune_snn) || !is.finite(prune_snn) || prune_snn < 0) {
    stop("'prune_snn' must be a finite non-negative value.", call. = FALSE)
  }

  k <- as.integer(k)
  if (prune_snn == 0) {
    return(0L)
  }
  threshold <- ((2 * k) * prune_snn) / (1 + prune_snn)
  max(0L, min(k, as.integer(ceiling(threshold - 1e-12))))
}

#' Construct Seurat-compatible SNN upper-triangle edges from ranked KNN.
#'
#' The returned edges are zero-based, row-major ordered, and unweighted for
#' construction of the undirected clustering graph.
#'
#' @param nn_idx_k Integer matrix of zero-based ranked KNN indices.
#'
#' @return An integer two-column matrix of upper-triangle edges.
#' @noRd
snn_upper_edge_matrix_from_knn_ranked_serial <- function(nn_idx_k) {
  if (!is.matrix(nn_idx_k) || !is.numeric(nn_idx_k) ||
      any(!is.finite(nn_idx_k)) ||
      any(nn_idx_k != floor(nn_idx_k))) {
    stop("'nn_idx_k' must be a finite integer matrix.", call. = FALSE)
  }

  n_cells <- nrow(nn_idx_k)
  k <- ncol(nn_idx_k)
  if (n_cells < 1L) {
    return(matrix(integer(), ncol = 2L))
  }
  if (any(nn_idx_k < 0L | nn_idx_k >= n_cells)) {
    stop("'nn_idx_k' contains an out-of-range neighbor index.", call. = FALSE)
  }

  min_shared <- min_shared_neighbors_for_prune(k)
  if (k == 0L) {
    return(matrix(integer(), ncol = 2L))
  }

  incidence <- Matrix::sparseMatrix(
    i = rep(seq_len(n_cells), each = k),
    j = as.integer(t(nn_idx_k)) + 1L,
    x = 1,
    dims = c(n_cells, n_cells)
  )
  upper <- Matrix::triu(Matrix::tcrossprod(incidence), k = 1L)
  triplets <- Matrix::summary(upper)
  if (NROW(triplets) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }
  triplets <- triplets[triplets$x >= min_shared, , drop = FALSE]
  if (NROW(triplets) == 0L) {
    return(matrix(integer(), ncol = 2L))
  }
  triplets <- triplets[order(triplets$i, triplets$j), , drop = FALSE]
  cbind(
    as.integer(triplets$i - 1L),
    as.integer(triplets$j - 1L)
  )
}


#' Construct Seurat-compatible SNN upper-triangle edges from ranked KNN.
#'
#' This dispatches independent row work to RcppParallel for non-trivial inputs.
#' @param nn_idx_k Integer matrix of zero-based ranked KNN indices.
#'
#' @return An integer two-column matrix of upper-triangle edges.
#' @noRd
snn_upper_edge_matrix_from_knn_ranked <- function(nn_idx_k) {
  if (!is.matrix(nn_idx_k) || !is.numeric(nn_idx_k) ||
      any(!is.finite(nn_idx_k)) ||
      any(nn_idx_k != floor(nn_idx_k))) {
    stop("'nn_idx_k' must be a finite integer matrix.", call. = FALSE)
  }

  n_cells <- nrow(nn_idx_k)
  k <- ncol(nn_idx_k)
  if (n_cells < 1L || k == 0L) {
    return(matrix(integer(), ncol = 2L))
  }
  if (any(nn_idx_k < 0L | nn_idx_k >= n_cells)) {
    stop("'nn_idx_k' contains an out-of-range neighbor index.", call. = FALSE)
  }

  min_shared <- min_shared_neighbors_for_prune(k)
  storage.mode(nn_idx_k) <- "integer"
  bcmp_snn_upper_edges_parallel_cpp(nn_idx_k, min_shared)
}
