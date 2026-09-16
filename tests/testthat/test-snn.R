test_that("SNN upper edges match the frozen Python fixture", {
  knn <- matrix(
    c(
      0L, 1L, 2L,
      1L, 0L, 3L,
      2L, 0L, 1L,
      3L, 1L, 0L,
      4L, 2L, 0L
    ),
    byrow = TRUE,
    ncol = 3L
  )

  expect_identical(
    bcmp:::snn_upper_edge_matrix_from_knn_ranked(knn),
    structure(
      c(
        0L, 0L, 0L, 0L, 1L, 1L, 1L, 2L, 2L, 3L,
        1L, 2L, 3L, 4L, 2L, 3L, 4L, 3L, 4L, 4L
      ),
      dim = c(10L, 2L)
    )
  )
})

test_that("SNN prune threshold retains the frozen Seurat semantics", {
  expect_identical(bcmp:::min_shared_neighbors_for_prune(0L), 0L)
  expect_identical(bcmp:::min_shared_neighbors_for_prune(3L), 1L)
  expect_identical(bcmp:::min_shared_neighbors_for_prune(15L), 2L)
  expect_identical(bcmp:::min_shared_neighbors_for_prune(15L, prune_snn = 0), 0L)
  expect_error(bcmp:::min_shared_neighbors_for_prune(3L, prune_snn = -1), "non-negative")
})

test_that("SNN pruning filters insufficient shared neighbors", {
  n_cells <- 30L
  k <- 15L
  knn <- t(vapply(
    0:(n_cells - 1L),
    function(cell) (cell + 0:(k - 1L)) %% n_cells,
    integer(k)
  ))
  edges <- bcmp:::snn_upper_edge_matrix_from_knn_ranked(knn)

  expect_true(any(edges[, 1L] == 0L & edges[, 2L] == 1L))
  expect_false(any(edges[, 1L] == 0L & edges[, 2L] == 15L))
})

test_that("SNN rejects invalid neighbor matrices", {
  expect_error(
    bcmp:::snn_upper_edge_matrix_from_knn_ranked(matrix(c(0, 1.5), nrow = 1L)),
    "integer matrix"
  )
  expect_error(
    bcmp:::snn_upper_edge_matrix_from_knn_ranked(matrix(c(0L, 2L), nrow = 1L)),
    "out-of-range"
  )
})


test_that("serial and parallel SNN preserve empty edge matrix contracts", {
  cases <- list(
    no_cells = matrix(integer(), nrow = 0L, ncol = 0L),
    no_neighbors = matrix(integer(), nrow = 3L, ncol = 0L),
    no_shared_neighbors = matrix(0:2, ncol = 1L)
  )
  expected <- matrix(integer(), nrow = 0L, ncol = 2L)
  for (name in names(cases)) {
    expect_identical(
      bcmp:::snn_upper_edge_matrix_from_knn_ranked_serial(cases[[name]]),
      expected, info = name
    )
    expect_identical(
      bcmp:::snn_upper_edge_matrix_from_knn_ranked(cases[[name]]),
      expected, info = name
    )
  }
})

test_that("serial and parallel SNN return no edges when all shared neighbors are pruned", {
  # Lines over a prime-sized grid share at most one point. Each row starts
  # with its own cell; eleven unique neighbors require two shared neighbors.
  side <- 11L
  n_cells <- side^2L
  knn <- t(vapply(0:(n_cells - 1L), function(cell) {
    slope <- cell %/% side
    intercept <- cell %% side - slope^2L
    x <- (slope + 0:(side - 1L)) %% side
    as.integer(x * side + (slope * x + intercept) %% side)
  }, integer(side)))
  incidence <- matrix(0L, nrow = n_cells, ncol = n_cells)
  incidence[cbind(rep(seq_len(n_cells), each = side), as.vector(t(knn)) + 1L)] <- 1L
  shared <- tcrossprod(incidence)
  expect_equal(range(shared[upper.tri(shared)]), c(0, 1))
  expect_identical(bcmp:::min_shared_neighbors_for_prune(side), 2L)

  expected <- matrix(integer(), nrow = 0L, ncol = 2L)
  expect_identical(bcmp:::snn_upper_edge_matrix_from_knn_ranked_serial(knn), expected)
  expect_identical(bcmp:::snn_upper_edge_matrix_from_knn_ranked(knn), expected)
})

test_that("parallel SNN matches the serial reference across thread counts", {
  set.seed(1)
  n_cells <- 192L
  k <- 15L
  knn <- t(vapply(
    seq_len(n_cells),
    function(cell) sample.int(n_cells, k) - 1L,
    integer(k)
  ))

  serial <- bcmp:::snn_upper_edge_matrix_from_knn_ranked_serial(knn)
  old_threads <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old_threads)) {
      Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
    } else {
      Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old_threads)
    }
  }, add = TRUE)

  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")
  one_thread <- bcmp:::bcmp_snn_upper_edges_parallel_cpp(
    knn,
    bcmp:::min_shared_neighbors_for_prune(k)
  )
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "4")
  four_threads <- bcmp:::bcmp_snn_upper_edges_parallel_cpp(
    knn,
    bcmp:::min_shared_neighbors_for_prune(k)
  )

  expect_identical(one_thread, serial)
  expect_identical(four_threads, serial)
  expect_identical(bcmp:::snn_upper_edge_matrix_from_knn_ranked(knn), serial)
})
