BCMP_NORMALIZATION_TARGET_SUM <- 1e4

validate_bcmp_count_matrix <- function(counts) {
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (!is.matrix(counts) && !inherits(counts, "Matrix")) {
    stop("`counts` must be a numeric matrix or Matrix sparse matrix.", call. = FALSE)
  }
  if (length(dim(counts)) != 2L || nrow(counts) < 1L || ncol(counts) < 1L) {
    stop("`counts` must have at least one cell and one feature.", call. = FALSE)
  }

  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.vector(counts)
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("`counts` must contain only finite count values.", call. = FALSE)
  }
  if (any(values < 0)) {
    stop("`counts` must contain non-negative count values.", call. = FALSE)
  }
  counts
}

normalize_log1p_counts <- function(counts, target_sum = BCMP_NORMALIZATION_TARGET_SUM) {
  counts <- validate_bcmp_count_matrix(counts)
  if (!is.numeric(target_sum) || length(target_sum) != 1L || is.na(target_sum) ||
      !is.finite(target_sum) || target_sum <= 0) {
    stop("`target_sum` must be one finite positive number.", call. = FALSE)
  }

  if (inherits(counts, "sparseMatrix")) {
    counts_sparse <- methods::as(counts, "generalMatrix")
    totals <- Matrix::rowSums(counts_sparse)
    factors <- ifelse(totals > 0, as.numeric(target_sum) / totals, 0)
    normalized <- Matrix::Diagonal(x = factors) %*% counts_sparse
    normalized@x <- log1p(normalized@x)
    return(normalized)
  }

  counts_dense <- as.matrix(counts)
  totals <- rowSums(counts_dense)
  factors <- ifelse(totals > 0, as.numeric(target_sum) / totals, 0)
  log1p(sweep(counts_dense, 1L, factors, "*"))
}

BCMP_BATCH_HVG_MIN_CELLS <- 3L

stabilize_duplicate_loess_predictor <- function(values, jitter_width = 1e-12) {
  values <- as.numeric(values)
  if (length(values) < 2L) {
    return(values)
  }
  order_idx <- order(values, method = "radix")
  stabilized <- values
  start <- 1L
  while (start <= length(values)) {
    end <- start
    while (end < length(values) && values[order_idx[end + 1L]] == values[order_idx[start]]) {
      end <- end + 1L
    }
    if (end > start) {
      stabilized[order_idx[start:end]] <- values[order_idx[start]] +
        seq(-jitter_width, jitter_width, length.out = end - start + 1L)
    }
    start <- end + 1L
  }
  stabilized
}

fit_bcmp_vst_loess <- function(log10_means, log10_variances, span = 0.3) {
  fit_once <- function(x) {
    stats::loess(
      log10_variances ~ x,
      span = span,
      degree = 2L,
      family = "gaussian"
    )
  }
  attempt <- tryCatch(fit_once(log10_means), error = identity)
  if (!inherits(attempt, "error")) {
    return(attempt)
  }
  stabilized <- stabilize_duplicate_loess_predictor(log10_means)
  if (identical(stabilized, as.numeric(log10_means))) {
    stop(conditionMessage(attempt), call. = FALSE)
  }
  retry <- tryCatch(fit_once(stabilized), error = identity)
  if (inherits(retry, "error")) {
    stop(conditionMessage(retry), call. = FALSE)
  }
  retry
}

as_bcmp_csc_count_matrix <- function(counts) {
  counts <- validate_bcmp_count_matrix(counts)
  if (!inherits(counts, "sparseMatrix")) {
    counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  }
  methods::as(methods::as(counts, "generalMatrix"), "CsparseMatrix")
}

select_hvgs_seurat_v5_vst <- function(raw_counts, var_names, n_top_genes, span = 0.3) {
  counts <- as_bcmp_csc_count_matrix(raw_counts)
  n_cells <- nrow(counts)
  n_genes <- ncol(counts)
  var_names <- as.character(var_names)
  n_top <- min(as.integer(n_top_genes), length(var_names))
  if (length(var_names) != n_genes || anyNA(var_names) || any(!nzchar(var_names)) || anyDuplicated(var_names)) {
    stop("`var_names` must contain one unique, non-empty name per count-matrix column.", call. = FALSE)
  }
  if (n_top < 1L) {
    stop("`n_top_genes` must be positive.", call. = FALSE)
  }
  if (n_cells < 2L) {
    stop("Seurat v5 VST requires at least 2 cells.", call. = FALSE)
  }

  sums <- Matrix::colSums(counts)
  means <- as.numeric(sums) / n_cells
  squared <- counts
  squared@x <- squared@x ^ 2
  sum_sq <- as.numeric(Matrix::colSums(squared))
  variances <- (sum_sq - (2 * means * sums) + (n_cells * means ^ 2)) / (n_cells - 1)
  variances <- pmax(variances, 0)

  expected_variance <- numeric(n_genes)
  non_constant <- variances > 0
  if (any(non_constant)) {
    fit <- fit_bcmp_vst_loess(
      log10(means[non_constant]),
      log10(variances[non_constant]),
      span = span
    )
    expected_variance[non_constant] <- 10 ^ fit$fitted
  }
  sd <- sqrt(expected_variance)
  vmax <- sqrt(n_cells)
  standardized_variance <- numeric(n_genes)
  for (gene_idx in seq_len(n_genes)) {
    if (sd[gene_idx] == 0 || !is.finite(sd[gene_idx])) {
      next
    }
    start <- counts@p[gene_idx] + 1L
    end <- counts@p[gene_idx + 1L]
    values <- if (start <= end) counts@x[start:end] else numeric()
    z <- (values - means[gene_idx]) / sd[gene_idx]
    z[z > vmax] <- vmax
    n_zero <- n_cells - length(values)
    zero_term <- ((0 - means[gene_idx]) / sd[gene_idx]) ^ 2
    standardized_variance[gene_idx] <- (sum(z ^ 2) + zero_term * n_zero) / (n_cells - 1)
  }

  candidate_idx <- which(means != 0)
  ranking <- candidate_idx[order(-standardized_variance[candidate_idx], method = "radix")]
  var_names[utils::head(ranking, n_top)]
}

seurat_v5_consensus_features <- function(features_by_layer, common_features, nfeatures = NULL) {
  if (!length(features_by_layer) || !length(common_features)) {
    return(character())
  }
  common_features <- as.character(common_features)
  positions <- list()
  for (features in features_by_layer) {
    features <- as.character(features)
    for (index in seq_along(features)) {
      feature <- features[[index]]
      if (feature %in% common_features) {
        positions[[feature]] <- c(positions[[feature]], index)
      }
    }
  }
  if (!length(positions)) {
    return(character())
  }
  features <- sort(names(positions))
  frequency <- vapply(features, function(feature) length(positions[[feature]]), integer(1))
  median_position <- vapply(features, function(feature) stats::median(positions[[feature]]), numeric(1))
  ranking <- order(-frequency, median_position, features, method = "radix")
  selected <- features[ranking]
  if (!is.null(nfeatures)) {
    selected <- utils::head(selected, as.integer(nfeatures))
  }
  selected
}

select_highly_variable_genes <- function(raw_counts, batch_labels, var_names, n_top_genes) {
  counts <- validate_bcmp_count_matrix(raw_counts)
  batch_labels <- as.character(batch_labels)
  if (length(batch_labels) != nrow(counts) || anyNA(batch_labels) || any(!nzchar(batch_labels))) {
    stop("`batch_labels` must contain one non-missing, non-empty label per cell.", call. = FALSE)
  }
  n_top <- min(as.integer(n_top_genes), ncol(counts))
  if (n_top < 1L) {
    stop("`n_top_genes` must be positive.", call. = FALSE)
  }
  batch_order <- unique(batch_labels)
  batch_sizes <- vapply(batch_order, function(batch) sum(batch_labels == batch), integer(1))
  small <- batch_sizes[batch_sizes < BCMP_BATCH_HVG_MIN_CELLS]
  if (length(small)) {
    preview <- paste(sprintf("%s=%d", names(small)[seq_len(min(5L, length(small)))], small[seq_len(min(5L, length(small)))]), collapse = ", ")
    suffix <- if (length(small) > 5L) ", ..." else ""
    stop(
      "Each batch must contain at least 3 cells for Seurat v5 split-layer HVG selection; small batches: ",
      preview,
      suffix,
      call. = FALSE
    )
  }
  features_by_layer <- lapply(batch_order, function(batch) {
    select_hvgs_seurat_v5_vst(
      raw_counts = counts[batch_labels == batch, , drop = FALSE],
      var_names = var_names,
      n_top_genes = n_top
    )
  })
  names(features_by_layer) <- paste0("counts.", batch_order)
  genes <- seurat_v5_consensus_features(features_by_layer, var_names, nfeatures = n_top)
  if (!length(genes)) {
    stop("No highly variable genes were selected by Seurat v5 split-layer consensus.", call. = FALSE)
  }
  genes
}

validate_bcmp_numeric_matrix <- function(matrix) {
  if (is.data.frame(matrix)) {
    matrix <- as.matrix(matrix)
  }
  if (!is.matrix(matrix) && !inherits(matrix, "Matrix")) {
    stop("`matrix` must be a numeric matrix or Matrix sparse matrix.", call. = FALSE)
  }
  values <- if (inherits(matrix, "sparseMatrix")) matrix@x else as.vector(matrix)
  if (!is.numeric(values) || any(!is.finite(values))) {
    stop("`matrix` must contain only finite numeric values.", call. = FALSE)
  }
  matrix
}

seurat_style_scale_dense <- function(matrix, scale_max = 10) {
  matrix <- validate_bcmp_numeric_matrix(matrix)
  if (!is.numeric(scale_max) || length(scale_max) != 1L || is.na(scale_max) ||
      !is.finite(scale_max) || scale_max <= 0) {
    stop("`scale_max` must be one finite positive number.", call. = FALSE)
  }
  dense <- as.matrix(matrix)
  means <- colMeans(dense)
  standard_deviations <- apply(dense, 2L, stats::sd)
  standard_deviations[!is.finite(standard_deviations) | standard_deviations == 0] <- 1
  scaled <- sweep(sweep(dense, 2L, means, "-"), 2L, standard_deviations, "/")
  scaled[scaled > scale_max] <- scale_max
  scaled[!is.finite(scaled)] <- 0
  scaled
}

run_bcmp_pca <- function(scaled_matrix, n_pcs) {
  scaled_matrix <- validate_bcmp_numeric_matrix(scaled_matrix)
  n_pcs <- as.integer(n_pcs)
  feasible_pcs <- min(nrow(scaled_matrix) - 1L, ncol(scaled_matrix))
  if (is.na(n_pcs) || n_pcs < 1L || feasible_pcs < n_pcs) {
    stop(
      "Requested partition_n_pcs but too few PCs are feasible for the current matrix.",
      call. = FALSE
    )
  }
  fit <- stats::prcomp(
    as.matrix(scaled_matrix),
    center = FALSE,
    scale. = FALSE,
    rank. = n_pcs,
    tol = 0
  )
  fit$x[, seq_len(n_pcs), drop = FALSE]
}

BCMP_UMAP_MAX_NEIGHBORS <- 30L
BCMP_UMAP_MAX_PCS <- 30L

run_bcmp_umap2d <- function(pca_embeddings, seed = 236L) {
  pca_embeddings <- validate_bcmp_numeric_matrix(pca_embeddings)
  n_cells <- nrow(pca_embeddings)
  n_pcs <- ncol(pca_embeddings)
  if (n_cells < 3L || n_pcs < 1L) {
    stop("BCMP debug UMAP requires at least 3 cells and 1 PCA dimension.", call. = FALSE)
  }
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop(
      "BCMP debug UMAP requires the suggested package `uwot`; install it with install.packages('uwot').",
      call. = FALSE
    )
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("`seed` must be one finite integer.", call. = FALSE)
  }

  if (n_pcs == 1L) {
    warning("Skipping BCMP debug UMAP because only one PCA dimension is available.",
            call. = FALSE)
    return(NULL)
  }

  dims_use <- seq_len(min(BCMP_UMAP_MAX_PCS, n_pcs))
  n_neighbors <- min(BCMP_UMAP_MAX_NEIGHBORS, n_cells - 1L)
  umap_embeddings <- with_bcmp_local_seed(seed, {
    uwot::umap(
      X = as.matrix(pca_embeddings[, dims_use, drop = FALSE]),
      n_threads = 1L,
      n_neighbors = as.integer(n_neighbors),
      n_components = 2L,
      metric = "cosine",
      n_epochs = NULL,
      learning_rate = 1,
      min_dist = 0.3,
      spread = 1,
      set_op_mix_ratio = 1,
      local_connectivity = 1L,
      repulsion_strength = 1,
      negative_sample_rate = 5L,
      a = NULL,
      b = NULL,
      fast_sgd = FALSE,
      approx_pow = FALSE,
      init = "spectral",
      verbose = FALSE,
      ret_model = FALSE
    )
  })
  colnames(umap_embeddings) <- c("UMAP_1", "UMAP_2")
  umap_embeddings
}

run_standard_bcmp_preprocess_from_counts <- function(
  raw_counts,
  batch_labels,
  var_names,
  partition_n_pcs = 30L,
  partition_n_hvg = 2000L,
  compute_umap = FALSE,
  seed = 236L
) {
  partition_n_pcs <- validate_bcmp_integer(partition_n_pcs, "partition_n_pcs", minimum = 1)
  partition_n_hvg <- validate_bcmp_integer(partition_n_hvg, "partition_n_hvg", minimum = 1)
  validate_bcmp_integer(seed, "seed")
  raw_counts <- validate_bcmp_count_matrix(raw_counts)
  if (nrow(raw_counts) < 3L) {
    stop("run_standard_bcmp_preprocess requires at least 3 cells.", call. = FALSE)
  }
  hvg_genes <- select_highly_variable_genes(
    raw_counts = raw_counts,
    batch_labels = batch_labels,
    var_names = var_names,
    n_top_genes = partition_n_hvg
  )
  hvg_idx <- match(hvg_genes, as.character(var_names))
  n_pcs <- as.integer(partition_n_pcs)
  if (length(hvg_idx) < n_pcs) {
    stop(
      "Requested partition_n_pcs=", n_pcs,
      " but only ", length(hvg_idx),
      " HVGs are available after feature selection.",
      call. = FALSE
    )
  }
  log_normalized <- normalize_log1p_counts(raw_counts)
  hvg_scaled <- seurat_style_scale_dense(log_normalized[, hvg_idx, drop = FALSE])
  pca_embeddings <- run_bcmp_pca(hvg_scaled, n_pcs)
  list(
    hvg_genes = hvg_genes,
    pca_embeddings = pca_embeddings,
    umap_embeddings = if (isTRUE(compute_umap)) run_bcmp_umap2d(pca_embeddings, seed = seed) else NULL
  )
}

run_bcmp_from_counts <- function(
  raw_counts,
  batch_labels,
  var_names,
  cell_ids = NULL,
  partition_n_pcs = 30L,
  partition_n_hvg = 2000L,
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
  validate_output_level(output_level)
  preprocess <- run_standard_bcmp_preprocess_from_counts(
    raw_counts = raw_counts,
    batch_labels = batch_labels,
    var_names = var_names,
    partition_n_pcs = partition_n_pcs,
    partition_n_hvg = partition_n_hvg,
    compute_umap = identical(output_level, "debug"),
    seed = seed
  )
  core <- bcmp_embedding(
    embedding = preprocess$pca_embeddings,
    batch_labels = batch_labels,
    cell_ids = cell_ids,
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
  list(core = core, preprocess = preprocess)
}
