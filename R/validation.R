#' Validate and canonicalize an embedding workflow input.
#'
#' The validator deliberately does not densify sparse inputs. It establishes
#' the input contract: cells are rows, batch labels align to
#' rows, and cell IDs are unique and reproducible.
#'
#' @param embedding A numeric matrix or supported sparse
#'   numeric matrix.
#' @param batch_labels One non-missing, non-empty batch label per cell.
#' @param cell_ids Optional explicit cell IDs.
#'
#' @return A list containing `embedding`, canonical `batch_labels`, and
#'   canonical `cell_ids`.
#' @noRd
validate_embedding_input <- function(embedding, batch_labels, cell_ids = NULL) {
  is_dense_matrix <- is.matrix(embedding) && is.numeric(embedding)
  is_sparse_matrix <- inherits(embedding, "sparseMatrix") &&
    "x" %in% methods::slotNames(embedding) && is.numeric(embedding@x)

  if (!(is_dense_matrix || is_sparse_matrix)) {
    stop(
      "`embedding` must be a numeric matrix or supported sparse numeric matrix.",
      call. = FALSE
    )
  }

  n_cells <- nrow(embedding)
  n_dimensions <- ncol(embedding)
  if (is.null(n_cells) || is.null(n_dimensions) || n_cells < 1L || n_dimensions < 1L) {
    stop("`embedding` must have at least one cell and one dimension.", call. = FALSE)
  }

  values <- if (is_sparse_matrix) embedding@x else unlist(embedding, use.names = FALSE)
  if (!all(is.finite(values))) {
    stop("`embedding` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.atomic(batch_labels) || length(batch_labels) != n_cells) {
    stop("`batch_labels` must contain one value for each embedding row.", call. = FALSE)
  }
  if (anyNA(batch_labels)) {
    stop("`batch_labels` must be non-missing, non-empty labels.", call. = FALSE)
  }
  batch_labels <- as.character(batch_labels)
  if (anyNA(batch_labels) || any(!nzchar(trimws(batch_labels)))) {
    stop("`batch_labels` must be non-missing, non-empty labels.", call. = FALSE)
  }

  if (is.null(cell_ids)) {
    cell_ids <- as.character(seq.int(0L, n_cells - 1L))
  } else {
    if (!is.atomic(cell_ids) || length(cell_ids) != n_cells) {
      stop("`cell_ids` must contain one value for each embedding row.", call. = FALSE)
    }
    if (anyNA(cell_ids)) {
      stop("`cell_ids` must be unique, non-missing, and non-empty.", call. = FALSE)
    }
    cell_ids <- as.character(cell_ids)
  }

  if (!has_valid_cell_ids(cell_ids, n_cells)) {
    stop("`cell_ids` must be unique, non-missing, and non-empty.", call. = FALSE)
  }

  list(
    embedding = embedding,
    batch_labels = batch_labels,
    cell_ids = cell_ids
  )
}

#' Check whether cell IDs satisfy the public embedding contract.
#'
#' @param cell_ids Candidate cell IDs.
#' @param n_cells Expected number of IDs.
#'
#' @return `TRUE` when the IDs are valid.
#' @noRd
has_valid_cell_ids <- function(cell_ids, n_cells) {
  !is.null(cell_ids) &&
    length(cell_ids) == n_cells &&
    !anyNA(cell_ids) &&
    all(nzchar(trimws(as.character(cell_ids)))) &&
    !anyDuplicated(cell_ids)
}

# Validate before coercion so fractional and out-of-range values are not truncated.
validate_bcmp_integer <- function(value, name, minimum = -.Machine$integer.max,
                                  allow_null = FALSE) {
  if (is.null(value) && allow_null) return(NULL)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < minimum || value > .Machine$integer.max ||
      value != trunc(value)) {
    stop("`", name, "` must be one finite integer in [", minimum, ", ",
         .Machine$integer.max, "]", if (allow_null) " or NULL" else "",
         ".", call. = FALSE)
  }
  as.integer(value)
}

validate_bcmp_partition_parameters <- function(
  k_min, k_max, min_batch_coverage, min_cells_in_domain_per_batch,
  max_underrepresentation_fold, seed
) {
  validate_bcmp_integer(k_min, "k_min", minimum = 1)
  validate_bcmp_integer(k_max, "k_max", minimum = 1, allow_null = TRUE)
  validate_bcmp_integer(min_batch_coverage, "min_batch_coverage",
                        minimum = 2, allow_null = TRUE)
  validate_bcmp_integer(min_cells_in_domain_per_batch,
                        "min_cells_in_domain_per_batch", minimum = 1)
  validate_bcmp_underrepresentation_fold(max_underrepresentation_fold)
  validate_bcmp_integer(seed, "seed")
  invisible(NULL)
}

validate_bcmp_underrepresentation_fold <- function(value) {
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
      is.na(value) || !is.finite(value) || value < 1) {
    stop("`max_underrepresentation_fold` must be one finite numeric value >= 1.",
         call. = FALSE)
  }
  as.numeric(value)
}
