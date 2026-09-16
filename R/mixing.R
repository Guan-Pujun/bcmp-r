BCMP_DEFAULT_MIN_BATCH_COVERAGE_FRAC <- 0.5
BCMP_DEFAULT_MIN_EFFECTIVE_CELLS_IN_DOMAIN_PER_BATCH <- 5L
BCMP_DEFAULT_MAX_UNDERREPRESENTATION_FOLD <- 10L
BCMP_DEFAULT_MAX_RESIDUAL_CELL_FRAC <- 0.05

required_batch_coverage <- function(n_batches, min_batch_coverage = NULL) {
  if (!is.numeric(n_batches) || length(n_batches) != 1L || is.na(n_batches) ||
      !is.finite(n_batches) || n_batches < 2L || n_batches != as.integer(n_batches)) {
    stop("Batch mixing requires at least two batches", call. = FALSE)
  }
  n_batches <- as.integer(n_batches)
  if (!is.null(min_batch_coverage)) {
    if (!is.numeric(min_batch_coverage) || length(min_batch_coverage) != 1L ||
        is.na(min_batch_coverage) || !is.finite(min_batch_coverage) ||
        min_batch_coverage != as.integer(min_batch_coverage) ||
        min_batch_coverage < 2L) {
      stop("min_batch_coverage must be at least 2 when supplied as a batch count", call. = FALSE)
    }
    if (min_batch_coverage > n_batches) {
      stop("min_batch_coverage cannot exceed the number of batches", call. = FALSE)
    }
    return(as.integer(min_batch_coverage))
  }
  as.integer(max(2L, ceiling(n_batches * BCMP_DEFAULT_MIN_BATCH_COVERAGE_FRAC)))
}

prepare_bcmp_batch_reference <- function(batch_labels) {
  if (!is.atomic(batch_labels) || !length(batch_labels)) {
    stop("prepare_bcmp_batch_reference requires at least one batch label", call. = FALSE)
  }
  batch_labels <- as.character(batch_labels)
  if (anyNA(batch_labels) || any(!nzchar(trimws(batch_labels)))) {
    stop("'batch_labels' must be non-missing, non-empty labels.", call. = FALSE)
  }
  batch_levels <- sort(unique(batch_labels))
  if (length(batch_levels) < 2L) {
    stop("Batch mixing requires at least two batches", call. = FALSE)
  }
  batch_idx <- match(batch_labels, batch_levels) - 1L
  batch_counts <- tabulate(batch_idx + 1L, nbins = length(batch_levels))
  list(
    batch_levels = batch_levels,
    batch_idx = as.integer(batch_idx),
    batch_counts = as.integer(batch_counts),
    p_global = as.numeric(batch_counts / length(batch_labels)),
    S = as.integer(length(batch_levels)),
    n = as.integer(length(batch_labels))
  )
}

evaluate_bcmp_batch_mixing_labels <- function(
  domain_labels,
  batch_labels = NULL,
  batch_reference = NULL,
  min_batch_coverage = NULL,
  min_cells_in_domain_per_batch = BCMP_DEFAULT_MIN_EFFECTIVE_CELLS_IN_DOMAIN_PER_BATCH,
  max_underrepresentation_fold = BCMP_DEFAULT_MAX_UNDERREPRESENTATION_FOLD,
  max_residual_cell_frac = BCMP_DEFAULT_MAX_RESIDUAL_CELL_FRAC
) {
  if (!is.atomic(domain_labels) || !length(domain_labels)) {
    stop("'domain_labels' must contain at least one label.", call. = FALSE)
  }
  domain_labels <- as.character(domain_labels)
  if (anyNA(domain_labels) || any(!nzchar(trimws(domain_labels)))) {
    stop("'domain_labels' must be non-missing, non-empty labels.", call. = FALSE)
  }
  if (is.null(batch_reference)) {
    if (is.null(batch_labels)) {
      stop("Either batch_labels or batch_reference must be supplied", call. = FALSE)
    }
    batch_reference <- prepare_bcmp_batch_reference(batch_labels)
  }
  if (length(domain_labels) != batch_reference$n) {
    stop("domain_labels length must match the encoded batch label length", call. = FALSE)
  }
  if (length(domain_labels) < 3L) {
    stop("Too few cells remaining for batch-mixing evaluation", call. = FALSE)
  }
  if (!is.numeric(min_cells_in_domain_per_batch) ||
      length(min_cells_in_domain_per_batch) != 1L ||
      is.na(min_cells_in_domain_per_batch) || !is.finite(min_cells_in_domain_per_batch) ||
      min_cells_in_domain_per_batch < 1L ||
      min_cells_in_domain_per_batch != as.integer(min_cells_in_domain_per_batch)) {
    stop("min_cells_in_domain_per_batch must be a positive integer", call. = FALSE)
  }
  max_underrepresentation_fold <- validate_bcmp_underrepresentation_fold(max_underrepresentation_fold)
  if (!is.numeric(max_residual_cell_frac) || length(max_residual_cell_frac) != 1L ||
      is.na(max_residual_cell_frac) || !is.finite(max_residual_cell_frac) ||
      max_residual_cell_frac < 0 || max_residual_cell_frac >= 1) {
    stop("max_residual_cell_frac must be a finite numeric fraction in [0, 1)", call. = FALSE)
  }

  domain_levels <- unique(domain_labels)
  domain_idx <- match(domain_labels, domain_levels) - 1L
  n_domains <- length(domain_levels)
  n_batches <- batch_reference$S
  required_coverage <- required_batch_coverage(n_batches, min_batch_coverage)
  counts <- matrix(
    tabulate(domain_idx + batch_reference$batch_idx * n_domains + 1L,
      nbins = n_domains * n_batches),
    nrow = n_domains,
    ncol = n_batches
  )
  domain_sizes <- rowSums(counts)
  first_seen <- match(domain_levels, domain_labels)
  domain_order <- order(-domain_sizes, first_seen)
  n_cells <- length(domain_labels)
  target_evaluated_cells <- min(
    max(ceiling((1 - max_residual_cell_frac) * n_cells), 1L),
    n_cells
  )
  evaluated <- rep(FALSE, n_domains)
  accumulated <- 0L
  for (domain_index in domain_order) {
    if (accumulated >= target_evaluated_cells) break
    evaluated[domain_index] <- TRUE
    accumulated <- accumulated + domain_sizes[domain_index]
  }

  thresholds <- outer(
    as.numeric(domain_sizes),
    batch_reference$p_global / max_underrepresentation_fold
  )
  effective <- counts >= as.integer(min_cells_in_domain_per_batch) & counts >= thresholds
  effective_count <- rowSums(effective)
  passes <- effective_count >= required_coverage
  domain_status <- ifelse(!evaluated, "residual", ifelse(passes, "pass", "fail"))
  evaluated_fail <- evaluated & !passes
  n_residual_cells <- sum(domain_sizes[!evaluated])

  per_domain <- data.frame(
    domain = domain_levels,
    domain_size = as.integer(domain_sizes),
    domain_size_frac = as.numeric(domain_sizes / n_cells),
    domain_first_seen_pos = as.integer(first_seen - 1L),
    effective_batch_count = as.integer(effective_count),
    domain_pass = domain_status,
    stringsAsFactors = FALSE
  )
  per_domain <- per_domain[domain_order, , drop = FALSE]
  mixing_status <- if (any(evaluated_fail)) "fail" else "pass"

  list(
    pass = identical(mixing_status, "pass"),
    mixing_status = mixing_status,
    fail_domain_count = as.integer(sum(evaluated_fail)),
    n_evaluated_domains = as.integer(sum(evaluated)),
    n_residual_domains = as.integer(sum(!evaluated)),
    n_residual_cells = as.integer(n_residual_cells),
    evaluated_cell_frac = as.numeric((n_cells - n_residual_cells) / n_cells),
    residual_cell_frac = as.numeric(n_residual_cells / n_cells),
    mixing_failure_reason = if (any(evaluated_fail)) "domain_mixing_fail" else "none",
    S = as.integer(n_batches),
    required_batch_coverage = as.integer(required_coverage),
    required_batch_coverage_source = if (is.null(min_batch_coverage)) "default_fraction" else "explicit",
    default_batch_coverage_fraction = BCMP_DEFAULT_MIN_BATCH_COVERAGE_FRAC,
    max_residual_cell_frac = as.numeric(max_residual_cell_frac),
    per_domain = per_domain
  )
}
