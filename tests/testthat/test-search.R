test_that("k-search records the Python-compatible low-k pass contract", {
  nn_idx <- matrix(
    c(0L, 1L, 2L,
      1L, 0L, 3L,
      2L, 0L, 1L,
      3L, 1L, 0L,
      4L, 2L, 0L),
    byrow = TRUE, ncol = 3L
  )
  result <- bcmp:::run_bcmp_k_search_from_knn(
    nn_idx,
    batch_labels = c("a", "b", "a", "b", "a"),
    cell_ids = as.character(0:4),
    k_low = 3L,
    k_high = 3L,
    min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0
  )

  expect_s3_class(result, "EmbeddingResult")
  expect_identical(result$selection$selected_k, 3L)
  expect_identical(result$selection$mixing_status, "pass")
  expect_identical(result$search_trace$phase, "k_low")
  expect_identical(names(result$labels), as.character(0:4))
})

test_that("debug k-search records candidates for evaluated scales", {
  nn_idx <- matrix(
    c(0L, 1L, 2L,
      1L, 0L, 3L,
      2L, 0L, 1L,
      3L, 1L, 0L,
      4L, 2L, 0L),
    byrow = TRUE, ncol = 3L
  )
  result <- bcmp:::run_bcmp_k_search_from_knn(
    nn_idx, c("a", "b", "a", "b", "a"), as.character(0:4),
    k_low = 3L, k_high = 3L, min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0, output_level = "debug"
  )

  expect_s3_class(result, "EmbeddingDebugResult")
  expect_identical(names(result$debug$candidates), "3")
  expect_identical(names(result$debug$candidates[[1L]]), c("partition", "domain_summary"))
})

test_that("bcmp_embedding runs the public embedding workflow", {
  embedding <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 2),
    c(3, 0),
    c(0, 5)
  )
  result <- bcmp_embedding(
    embedding,
    batch_labels = c("a", "b", "a", "b", "a"),
    k_min = 3L,
    k_max = 3L,
    min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0
  )

  expect_s3_class(result, "EmbeddingResult")
  expect_identical(result$selection$selected_k, 3L)
  expect_identical(names(result$labels), as.character(0:4))
})


test_that("search returns the upper-bound diagnostic when neither bound passes", {
  # Two separated, batch-pure groups remain unmixed throughout this range.
  embedding <- cbind(c((0:9)^2 / 10, 100 + (0:9)^2 / 9), 0)
  result <- bcmp_embedding(
    embedding, rep(c("a", "b"), each = 10L),
    k_min = 3L, k_max = 5L, min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0, output_level = "debug", verbose = FALSE
  )

  expect_identical(result$search_trace$k, c(3L, 5L))
  expect_identical(result$search_trace$phase, c("k_low", "k_high"))
  expect_identical(result$search_trace$mixing_status, c("fail", "fail"))
  expect_identical(result$debug$search$status, "diagnostic_no_feasible_k")
  expect_identical(result$selection$selected_k, 5L)
  expect_identical(result$selection$mixing_status, "fail")
  expect_identical(result$labels, result$debug$candidates[["5"]]$partition)
})

test_that("search advances past failing midpoints and keeps the best passing scale", {
  # Larger neighborhoods connect the same two groups and permit mixing.
  embedding <- cbind(c((0:9)^2 / 10, 100 + (0:9)^2 / 9), 0)
  result <- bcmp_embedding(
    embedding, rep(c("a", "b"), each = 10L),
    k_min = 3L, k_max = 19L, min_cells_in_domain_per_batch = 1L,
    max_residual_cell_frac = 0, output_level = "debug", verbose = FALSE
  )

  expect_identical(result$search_trace$k, c(3L, 19L, 11L, 7L, 9L, 10L))
  expect_identical(
    result$search_trace$phase, c("k_low", "k_high", rep("k_mid", 4L))
  )
  expect_identical(
    result$search_trace$mixing_status,
    c("fail", "pass", "pass", "fail", "fail", "fail")
  )
  expect_identical(result$debug$search$status, "binary_search_best_pass")
  expect_identical(result$selection$selected_k, 11L)
  expect_identical(result$selection$mixing_status, "pass")
  expect_identical(result$labels, result$debug$candidates[["11"]]$partition)
})

test_that("verbose emits Python-style k-search records", {
  nn_idx <- matrix(
    c(0L, 1L, 2L,
      1L, 0L, 3L,
      2L, 0L, 1L,
      3L, 1L, 0L,
      4L, 2L, 0L),
    byrow = TRUE, ncol = 3L
  )
  call <- function(verbose) {
    bcmp:::run_bcmp_k_search_from_knn(
      nn_idx, c("a", "b", "a", "b", "a"), as.character(0:4),
      k_low = 3L, k_high = 3L, min_cells_in_domain_per_batch = 1L,
      max_residual_cell_frac = 0, verbose = verbose
    )
  }

  expect_message(call(TRUE), "\\[iter 1\\]\\[k_low\\] k=3")
  expect_silent(call(FALSE))
})
