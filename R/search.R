format_bcmp_progress_line <- function(record) {
  sprintf(
    "[iter %d][%s] k=%d | mixing_status=%s | n_domains=%d | n_failed_domains=%d | n_residual_domains=%d | actual_residual_cell_frac=%.4f",
    as.integer(record$iter_id[[1L]]),
    as.character(record$phase[[1L]]),
    as.integer(record$k[[1L]]),
    as.character(record$mixing_status[[1L]]),
    as.integer(record$n_domains[[1L]]),
    as.integer(record$n_failed_domains[[1L]]),
    as.integer(record$n_residual_domains[[1L]]),
    as.numeric(record$actual_residual_cell_frac[[1L]])
  )
}

emit_bcmp_progress_line <- function(record, verbose) {
  if (isTRUE(verbose)) {
    message(format_bcmp_progress_line(record))
  }
  invisible(NULL)
}

format_bcmp_selected_progress <- function(record) {
  sprintf(
    "[selected] k=%d | mixing_status=%s | n_domains=%d | n_failed_domains=%d | n_residual_domains=%d | actual_residual_cell_frac=%.4f",
    as.integer(record$k[[1L]]),
    as.character(record$mixing_status[[1L]]),
    as.integer(record$n_domains[[1L]]),
    as.integer(record$n_failed_domains[[1L]]),
    as.integer(record$n_residual_domains[[1L]]),
    as.numeric(record$actual_residual_cell_frac[[1L]])
  )
}

labels_at_k <- function(nn_idx_max, k, seed = 236L) {
  if (!is.matrix(nn_idx_max) || k < 1L || k > ncol(nn_idx_max)) {
    stop("Invalid KNN slice for requested k.", call. = FALSE)
  }
  edges <- snn_upper_edge_matrix_from_knn_ranked(nn_idx_max[, seq_len(k), drop = FALSE])
  cluster_labels_from_edge_matrix(nrow(nn_idx_max), edges, seed = seed)
}

new_bcmp_trace <- function(rows) {
  out <- data.frame(
    iter_id = integer(),
    phase = character(),
    k = integer(),
    mixing_status = character(),
    n_domains = integer(),
    n_failed_domains = integer(),
    n_residual_domains = integer(),
    actual_residual_cell_frac = numeric(),
    stringsAsFactors = FALSE
  )
  if (length(rows)) {
    out <- do.call(rbind, rows)
  }
  out[, BCMP_SEARCH_TRACE_FIELDS, drop = FALSE]
}

run_bcmp_k_search_from_knn <- function(
  nn_idx_max,
  batch_labels,
  cell_ids,
  k_low,
  k_high,
  min_batch_coverage = NULL,
  min_cells_in_domain_per_batch = BCMP_DEFAULT_MIN_EFFECTIVE_CELLS_IN_DOMAIN_PER_BATCH,
  max_underrepresentation_fold = BCMP_DEFAULT_MAX_UNDERREPRESENTATION_FOLD,
  max_residual_cell_frac = BCMP_DEFAULT_MAX_RESIDUAL_CELL_FRAC,
  seed = 236L,
  output_level = "standard",
  verbose = TRUE
) {
  if (!is.matrix(nn_idx_max) || nrow(nn_idx_max) != length(batch_labels) ||
      nrow(nn_idx_max) != length(cell_ids) || k_low > k_high ||
      k_low < 1L || k_high > ncol(nn_idx_max)) {
    stop("Invalid KNN search inputs.", call. = FALSE)
  }
  batch_reference <- prepare_bcmp_batch_reference(batch_labels)
  trace_rows <- list()
  debug_candidates <- list()
  iter_id <- 0L

  evaluate_k <- function(k, phase) {
    iter_id <<- iter_id + 1L
    labels_num <- labels_at_k(nn_idx_max, k, seed = seed)
    labels <- as.character(labels_num)
    mix <- evaluate_bcmp_batch_mixing_labels(
      labels,
      batch_reference = batch_reference,
      min_batch_coverage = min_batch_coverage,
      min_cells_in_domain_per_batch = min_cells_in_domain_per_batch,
      max_underrepresentation_fold = max_underrepresentation_fold,
      max_residual_cell_frac = max_residual_cell_frac
    )
    trace_row <- data.frame(
      iter_id = as.integer(iter_id),
      phase = phase,
      k = as.integer(k),
      mixing_status = mix$mixing_status,
      n_domains = as.integer(nrow(mix$per_domain)),
      n_failed_domains = mix$fail_domain_count,
      n_residual_domains = mix$n_residual_domains,
      actual_residual_cell_frac = mix$residual_cell_frac,
      stringsAsFactors = FALSE
    )
    trace_rows[[length(trace_rows) + 1L]] <<- trace_row
    emit_bcmp_progress_line(trace_row, verbose)
    if (identical(output_level, "debug")) {
      summary <- mix$per_domain[
        c("domain", "domain_size", "domain_size_frac", "effective_batch_count", "domain_pass")
      ]
      names(summary)[names(summary) == "domain_size_frac"] <- "domain_cell_frac"
      names(summary)[names(summary) == "domain_pass"] <- "domain_mixing_status"
      debug_candidates[[as.character(k)]] <<- list(
        partition = stats::setNames(labels, cell_ids),
        domain_summary = summary
      )
    }
    list(k = as.integer(k), labels = labels, mix = mix)
  }

  low <- evaluate_k(k_low, "k_low")
  if (k_low == k_high || low$mix$pass) {
    chosen <- low
    search_status <- if (low$mix$pass) "k_low_pass" else "diagnostic_no_feasible_k"
  } else {
    high <- evaluate_k(k_high, "k_high")
    if (!high$mix$pass) {
      chosen <- high
      search_status <- "diagnostic_no_feasible_k"
    } else {
      best <- high
      low_k <- k_low
      high_k <- k_high
      while (high_k - low_k > 1L) {
        mid_k <- as.integer(floor((low_k + high_k) / 2))
        mid <- evaluate_k(mid_k, "k_mid")
        if (mid$mix$pass) {
          best <- mid
          high_k <- mid_k
        } else {
          low_k <- mid_k
        }
      }
      chosen <- best
      search_status <- "binary_search_best_pass"
    }
  }

  selected_record <- trace_rows[[which(vapply(
    trace_rows,
    function(record) as.integer(record$k[[1L]]) == chosen$k,
    logical(1L)
  ))[[1L]]]]
  if (isTRUE(verbose)) {
    message(format_bcmp_selected_progress(selected_record))
  }

  selection <- new_bcmp_selection(
    selected_k = chosen$k,
    n_domains = as.integer(nrow(chosen$mix$per_domain)),
    mixing_status = chosen$mix$mixing_status,
    n_evaluated_domains = chosen$mix$n_evaluated_domains,
    n_residual_domains = chosen$mix$n_residual_domains,
    n_residual_cells = chosen$mix$n_residual_cells,
    residual_cell_frac = chosen$mix$residual_cell_frac
  )
  trace <- new_bcmp_trace(trace_rows)
  labels <- stats::setNames(chosen$labels, cell_ids)
  if (identical(output_level, "debug")) {
    debug <- list(
      parameters = list(
        min_batch_coverage = min_batch_coverage,
        min_cells_in_domain_per_batch = as.integer(min_cells_in_domain_per_batch),
        max_underrepresentation_fold = as.numeric(max_underrepresentation_fold),
        max_residual_cell_frac = max_residual_cell_frac,
        seed = as.integer(seed)
      ),
      search = list(
        k_low = as.integer(k_low),
        k_high = as.integer(k_high),
        status = search_status,
        required_batch_coverage = chosen$mix$required_batch_coverage,
        required_batch_coverage_source = chosen$mix$required_batch_coverage_source
      ),
      candidates = debug_candidates
    )
    return(new_embedding_debug_result(labels, selection, trace, debug))
  }
  new_embedding_result(labels, selection, trace)
}

run_bcmp_partition_search <- function(
  embedding,
  batch_labels,
  cell_ids,
  k_min = 3L,
  k_max = NULL,
  batch_frac_threshold_for_k_max = 0.01,
  min_batch_coverage = NULL,
  min_cells_in_domain_per_batch = BCMP_DEFAULT_MIN_EFFECTIVE_CELLS_IN_DOMAIN_PER_BATCH,
  max_underrepresentation_fold = BCMP_DEFAULT_MAX_UNDERREPRESENTATION_FOLD,
  max_residual_cell_frac = BCMP_DEFAULT_MAX_RESIDUAL_CELL_FRAC,
  seed = 236L,
  output_level = "standard",
  verbose = TRUE
) {
  validate_bcmp_partition_parameters(
    k_min, k_max, min_batch_coverage, min_cells_in_domain_per_batch,
    max_underrepresentation_fold, seed
  )
  bounds <- derive_k_search_bounds(
    batch_labels, nrow(embedding), k_min, k_max, batch_frac_threshold_for_k_max
  )
  if (bounds$k_low > bounds$k_high) {
    stop("Invalid k range", call. = FALSE)
  }
  nn_idx_max <- get_knn_neighbor_once(embedding, bounds$k_high, seed)
  run_bcmp_k_search_from_knn(
    nn_idx_max = nn_idx_max,
    batch_labels = batch_labels,
    cell_ids = cell_ids,
    k_low = bounds$k_low,
    k_high = bounds$k_high,
    min_batch_coverage = min_batch_coverage,
    min_cells_in_domain_per_batch = min_cells_in_domain_per_batch,
    max_underrepresentation_fold = max_underrepresentation_fold,
    max_residual_cell_frac = max_residual_cell_frac,
    seed = seed,
    output_level = output_level,
    verbose = verbose
  )
}
