test_that("a sole positive domain candidate always receives the singleton", {
  edges <- rbind(c(0L, 1L), c(1L, 2L))
  for (candidate in c(1L, 2L, 7L)) {
    for (seed in c(1:10, 236L)) {
      expect_identical(
        bcmp:::group_singletons_from_edges(
          c(0L, candidate, candidate), 3L, edges, seed = seed
        ),
        c(0L, 0L, 0L)
      )
    }
  }
})

test_that("singleton ties select only eligible domains and are reproducible", {
  labels <- c(0L, 2L, 2L, 7L, 7L, 9L, 9L)
  edges <- rbind(c(0L, 1L), c(0L, 3L), c(1L, 2L),
                 c(3L, 4L), c(5L, 6L))
  for (seed in c(1:10, 236L)) {
    result <- bcmp:::group_singletons_from_edges(labels, 7L, edges, seed)
    expect_true(result[1L] %in% result[c(2L, 4L)])
    expect_false(result[1L] == result[6L])
    expect_identical(result[-1L], c(0L, 0L, 1L, 1L, 2L, 2L))
    expect_identical(
      result, bcmp:::group_singletons_from_edges(labels, 7L, edges, seed)
    )
  }
})

test_that("singleton grouping follows the Python graph-connectivity rule", {
  edges <- rbind(
    c(0L, 1L),
    c(2L, 3L),
    c(0L, 4L)
  )
  expect_identical(
    bcmp:::group_singletons_from_edges(
      labels = c(0L, 0L, 1L, 1L, 2L),
      n_cells = 5L,
      edges = edges,
      seed = 236L
    ),
    c(0L, 0L, 1L, 1L, 0L)
  )
})

test_that("singleton grouping restores the caller RNG state", {
  set.seed(17L)
  first <- runif(1L)
  edges <- rbind(c(0L, 1L), c(2L, 3L), c(0L, 4L))
  bcmp:::group_singletons_from_edges(
    labels = c(0L, 0L, 1L, 1L, 2L),
    n_cells = 5L,
    edges = edges,
    seed = 236L
  )
  after <- runif(1L)

  set.seed(17L)
  expected_first <- runif(1L)
  expected_after <- runif(1L)
  expect_identical(first, expected_first)
  expect_identical(after, expected_after)
})

test_that("all-singleton initial clustering fails explicitly", {
  expect_error(
    bcmp:::group_singletons_from_edges(
      labels = c(0L, 1L, 2L),
      n_cells = 3L,
      edges = matrix(integer(), ncol = 2L)
    ),
    "All clusters are singletons"
  )
})

test_that("Leiden preserves a disconnected-clique partition", {
  edges <- rbind(
    c(0L, 1L), c(0L, 2L), c(1L, 2L),
    c(3L, 4L), c(3L, 5L), c(4L, 5L)
  )
  labels <- bcmp:::cluster_labels_from_edge_matrix(6L, edges, seed = 236L)

  expect_true(all(labels[1:3] == labels[1L]))
  expect_true(all(labels[4:6] == labels[4L]))
  expect_false(labels[1L] == labels[4L])
})

test_that("cluster core validates upper-triangle edges", {
  expect_error(
    bcmp:::cluster_labels_from_edge_matrix(3L, matrix(c(1L, 0L), ncol = 2L)),
    "upper-triangle"
  )
})
