test_that("matrix input receives Python-compatible default cell IDs", {
  embedding <- matrix(c(1, 2, 3, 4), nrow = 2)

  validated <- bcmp:::validate_embedding_input(
    embedding = embedding,
    batch_labels = c("batch_a", "batch_b")
  )

  expect_identical(validated$cell_ids, c("0", "1"))
  expect_identical(validated$batch_labels, c("batch_a", "batch_b"))
  expect_identical(validated$embedding, embedding)
})

test_that("data frames are rejected at the public embedding entry point", {
  x <- data.frame(pc_1 = 1:20, pc_2 = (1:20)^2)
  expect_error(bcmp_embedding(x, rep(c("a", "b"), each = 10)), "numeric matrix")
})

test_that("invalid embedding and alignment inputs fail explicitly", {
  expect_error(
    bcmp:::validate_embedding_input(data.frame(value = c("a", "b")), c("a", "b")),
    "numeric matrix"
  )
  expect_error(
    bcmp:::validate_embedding_input(matrix(c(1, Inf), nrow = 2), c("a", "b")),
    "finite"
  )
  expect_error(
    bcmp:::validate_embedding_input(matrix(1:4, nrow = 2), "batch_a"),
    "one value"
  )
  expect_error(
    bcmp:::validate_embedding_input(
      matrix(1:4, nrow = 2),
      c("batch_a", "batch_b"),
      cell_ids = c("same", "same")
    ),
    "unique"
  )
})

test_that("public embedding input rejects empty dimensions", {
  expect_error(
    bcmp_embedding(matrix(numeric(), 0L, 2L), character()),
    "at least one cell and one dimension"
  )
  expect_error(
    bcmp_embedding(matrix(numeric(), 4L, 0L), c("a", "a", "b", "b")),
    "at least one cell and one dimension"
  )
})

test_that("public embedding input rejects missing metadata and misaligned cell IDs", {
  embedding <- matrix(seq_len(8L), nrow = 4L)
  for (label in c(NA_character_, "", " ")) {
    expect_error(
      bcmp_embedding(embedding, c("a", "a", "b", label)),
      "non-missing, non-empty labels"
    )
  }
  expect_error(
    bcmp_embedding(embedding, c("a", "a", "b", "b"), cell_ids = c("c1", "c2", "c3")),
    "one value for each embedding row"
  )
})

test_that("numeric missing metadata is rejected before character conversion", {
  embedding <- matrix(seq_len(8L), nrow = 4L)
  for (missing in list(NA_real_, NaN)) {
    expect_error(
      bcmp_embedding(embedding, c(1, 1, 2, missing)),
      "non-missing, non-empty labels"
    )
    expect_error(
      bcmp_embedding(embedding, c(1, 1, 2, 2), cell_ids = c(1, 2, 3, missing)),
      "unique, non-missing, and non-empty"
    )
  }
})

test_that("numeric labels and literal NaN strings remain valid metadata", {
  embedding <- matrix(seq_len(8L), nrow = 4L)
  for (values in list(seq_len(4L), c("a", "b", "c", "NaN"))) {
    validated <- bcmp:::validate_embedding_input(embedding, values, values)
    expect_identical(validated$batch_labels, as.character(values))
    expect_identical(validated$cell_ids, as.character(values))
    expect_identical(
      unname(bcmp:::validate_bcmp_object_batches(values, as.character(1:4), "batch")),
      as.character(values)
    )
  }
})

test_that("only documented output levels are accepted", {
  expect_silent(bcmp:::validate_output_level("standard"))
  expect_silent(bcmp:::validate_output_level("debug"))
  expect_error(bcmp:::validate_output_level("minimal"), "standard")
})

test_that("result constructors retain Python-compatible classes and field order", {
  selection <- bcmp:::new_bcmp_selection(
    selected_k = 3L,
    n_domains = 2L,
    mixing_status = "pass",
    n_evaluated_domains = 2L,
    n_residual_domains = 0L,
    n_residual_cells = 0L,
    residual_cell_frac = 0
  )
  search_trace <- data.frame(
    iter_id = 1L,
    phase = "k_low",
    k = 3L,
    mixing_status = "pass",
    n_domains = 2L,
    n_failed_domains = 0L,
    n_residual_domains = 0L,
    actual_residual_cell_frac = 0
  )

  result <- bcmp:::new_embedding_debug_result(
    labels = c(cell_a = "Domain_1", cell_b = "Domain_2"),
    selection = selection,
    search_trace = search_trace,
    debug = list()
  )

  expect_identical(names(selection), bcmp:::BCMP_SELECTION_FIELDS)
  expect_identical(names(search_trace), bcmp:::BCMP_SEARCH_TRACE_FIELDS)
  expect_identical(names(result), c("labels", "selection", "search_trace", "debug"))
  expect_identical(class(result), c("EmbeddingDebugResult", "EmbeddingResult"))
})

test_that("k-search bounds match the Python-v1 eligible batch rule", {
  expect_identical(
    bcmp:::derive_k_search_bounds(
      batch_labels = c(rep("a", 6), rep("b", 4), rep("c", 5)),
      n_total = 15L,
      k_min = 3L
    ),
    list(S = 3L, k_low = 3L, k_high = 3L)
  )
  expect_identical(
    bcmp:::derive_k_search_bounds(
      batch_labels = c("tiny", rep("a", 6), rep("b", 5)),
      n_total = 12L,
      k_min = 3L,
      batch_frac_threshold_for_k_max = 0.1
    ),
    list(S = 3L, k_low = 3L, k_high = 4L)
  )
})

test_that("k-search bounds use a strict threshold and cap explicit k_max", {
  expect_error(
    bcmp:::derive_k_search_bounds(
      batch_labels = c(rep("a", 3), rep("b", 3)),
      n_total = 6L,
      batch_frac_threshold_for_k_max = 0.5
    ),
    "no batches exceed 50.00%"
  )
  expect_identical(
    bcmp:::derive_k_search_bounds(
      batch_labels = c(rep("a", 8), rep("b", 8)),
      n_total = 16L,
      k_max = 99L
    ),
    list(S = 2L, k_low = 3L, k_high = 15L)
  )
  expect_error(
    bcmp:::derive_k_search_bounds(c("a", "b"), 2L, batch_frac_threshold_for_k_max = 1),
    "fraction"
  )
})

test_that("integer scalar validation handles range, type and missing values", {
  expect_identical(bcmp:::validate_bcmp_integer(3, "k_min", minimum = 1), 3L)
  expect_null(bcmp:::validate_bcmp_integer(NULL, "k_max", allow_null = TRUE))
  for (x in list(2.5, NA_real_, Inf, NaN, numeric(), c(2, 3), "3", TRUE,
                 .Machine$integer.max + 1)) {
    expect_error(bcmp:::validate_bcmp_integer(x, "k_min", minimum = 1), "k_min")
  }
  expect_error(bcmp:::validate_bcmp_integer(0, "k_min", minimum = 1), "k_min")
})

test_that("embedding public integer arguments reject fractional values", {
  args <- list(embedding = matrix(seq_len(60), nrow = 20),
               batch_labels = rep(c("a", "b"), each = 10),
               k_max = 3L, verbose = FALSE)
  for (name in c("k_min", "k_max", "min_batch_coverage",
                 "min_cells_in_domain_per_batch", "seed")) {
    invalid <- args
    invalid[[name]] <- 2.5
    expect_error(do.call(bcmp_embedding, invalid), name)
  }
})

test_that("integer preprocessing settings reject fractional values", {
  for (name in c("partition_n_pcs", "partition_n_hvg")) {
    args <- list(raw_counts = matrix(1:120, nrow = 20),
                 batch_labels = rep(c("a", "b"), each = 10),
                 var_names = paste0("g", 1:6), partition_n_pcs = 2L,
                 partition_n_hvg = 4L)
    args[[name]] <- 2.5
    expect_error(do.call(bcmp:::run_standard_bcmp_preprocess_from_counts, args), name)
  }
})

test_that("fold accepts finite real scalars >= 1 and preserves debug precision", {
  args <- list(embedding = matrix(seq_len(120), nrow = 40),
               batch_labels = rep(c("a", "b"), each = 20),
               k_min = 3L, k_max = 3L, output_level = "debug", verbose = FALSE)
  for (fold in c(1, 2.5, 3, 3.33, 10)) {
    result <- do.call(bcmp_embedding, c(args, list(max_underrepresentation_fold = fold)))
    expect_identical(result$debug$parameters$max_underrepresentation_fold, fold)
  }
  for (fold in list(0, -1, 0.5, NA_real_, NaN, Inf, -Inf, "2.5", TRUE,
                    2+1i, numeric(), c(2, 3), NULL)) {
    expect_error(do.call(bcmp_embedding, c(args, list(max_underrepresentation_fold = fold))),
                 "max_underrepresentation_fold")
  }
})

test_that("explicit k_max bypasses automatic batch eligibility", {
  batches <- rep(c("a", "b"), each = 20)
  for (threshold in c(0, 0.5, 0.9)) {
    expect_identical(
      bcmp:::derive_k_search_bounds(batches, 40L, k_max = 3L,
        batch_frac_threshold_for_k_max = threshold),
      list(S = 2L, k_low = 3L, k_high = 3L))
  }
  expect_identical(
    bcmp:::derive_k_search_bounds(batches, 40L, k_max = 99L,
      batch_frac_threshold_for_k_max = 0.5),
    list(S = 2L, k_low = 3L, k_high = 39L))
  expect_error(
    bcmp:::derive_k_search_bounds(batches, 40L,
      batch_frac_threshold_for_k_max = 0.5), "no batches exceed")
  expect_error(
    bcmp:::derive_k_search_bounds(batches, 40L, k_max = 3L,
      batch_frac_threshold_for_k_max = 1), "fraction")
  x <- matrix(seq_len(120), nrow = 40)
  baseline <- bcmp_embedding(x, batches, k_max = 3L, verbose = FALSE)
  explicit <- bcmp_embedding(x, batches, k_max = 3L,
    batch_frac_threshold_for_k_max = 0.5, verbose = FALSE)
  expect_identical(explicit, baseline)
})
