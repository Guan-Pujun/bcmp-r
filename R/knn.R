#' Build one Annoy KNN index for a fixed embedding.
#'
#' Uses Euclidean distance, 50 trees, and a local seed. The search setting -1
#' uses Annoy's default search budget (requested neighbors times tree count),
#' not an exhaustive search.
#'
#' @param embeddings Numeric cells-by-dimensions embedding.
#' @param k_max Requested number of neighbors.
#' @param seed Local Annoy seed.
#'
#' @return An integer matrix of zero-based neighbor indices.
#' @noRd
get_knn_neighbor_once <- function(embeddings, k_max, seed = 236L) {
  is_dense_matrix <- is.matrix(embeddings) && is.numeric(embeddings)
  is_sparse_matrix <- inherits(embeddings, "sparseMatrix") &&
    "x" %in% methods::slotNames(embeddings) && is.numeric(embeddings@x)
  if (!(is_dense_matrix || is_sparse_matrix)) {
    stop("'embeddings' must be a finite numeric matrix or supported sparse numeric matrix.", call. = FALSE)
  }
  values <- if (is_sparse_matrix) embeddings@x else embeddings
  if (any(!is.finite(values))) {
    stop("'embeddings' must contain only finite values.", call. = FALSE)
  }

  n_cells <- nrow(embeddings)
  if (n_cells < 3L) {
    stop("At least 3 cells are required for KNN search", call. = FALSE)
  }
  if (ncol(embeddings) < 1L) {
    stop("'embeddings' must contain at least one dimension.", call. = FALSE)
  }
  # Annoy uses float32, including squared distances and up to 200 weighted
  # centroid updates. Bound 4 * dimensions * magnitude^2 * 201 with 4x headroom.
  float32_max <- (2 - 2^-23) * 2^127
  max_magnitude <- sqrt(float32_max / (16 * 201 * ncol(embeddings)))
  if (any(abs(values) > max_magnitude)) {
    stop("'embeddings' magnitude is too large for safe float32 Annoy computation.",
         call. = FALSE)
  }
  if (!is.numeric(k_max) || length(k_max) != 1L || is.na(k_max) ||
      !is.finite(k_max) || k_max != as.integer(k_max)) {
    stop("'k_max' must be one finite integer.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("'seed' must be one finite integer.", call. = FALSE)
  }

  k_use <- min(as.integer(k_max), n_cells - 1L)
  if (k_use < 1L) {
    stop("k_max must be at least 1", call. = FALSE)
  }

  index <- methods::new(RcppAnnoy::AnnoyEuclidean, ncol(embeddings))
  index$setSeed(as.integer(seed))
  for (cell_index in seq_len(n_cells)) {
    index$addItem(cell_index - 1L, as.numeric(embeddings[cell_index, ]))
  }
  index$build(50L)

  out <- matrix(NA_integer_, nrow = n_cells, ncol = k_use)
  for (cell_index in seq_len(n_cells)) {
    neighbors <- index$getNNsByVectorList(
      as.numeric(embeddings[cell_index, ]),
      k_use,
      -1L,
      FALSE
    )$item
    if (length(neighbors) != k_use) {
      stop(
        sprintf(
          "Annoy returned %d neighbors for cell %d, expected %d",
          length(neighbors),
          cell_index - 1L,
          k_use
        ),
        call. = FALSE
      )
    }
    out[cell_index, ] <- as.integer(neighbors)
  }
  out
}
