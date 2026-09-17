test_that("fixed KNN matches the frozen Python Annoy fixture", {
  embedding <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 2),
    c(3, 0),
    c(0, 5)
  )

  neighbors <- bcmp:::get_knn_neighbor_once(embedding, k_max = 3L, seed = 236L)
  expect_identical(dim(neighbors), c(5L, 3L))
  expect_type(neighbors, "integer")
  expect_identical(neighbors[, 1L], 0:4)
  expect_true(all(apply(neighbors, 1L, anyDuplicated) == 0L))

  expect_identical(
    neighbors,
    structure(
      c(0L, 1L, 2L, 3L, 4L, 1L, 0L, 0L, 1L, 2L, 2L, 3L, 1L, 0L, 0L),
      dim = c(5L, 3L)
    )
  )
})

test_that("fixed KNN handles equidistant neighbors without requiring tie order", {
  embedding <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  neighbors <- bcmp:::get_knn_neighbor_once(embedding, k_max = 3L, seed = 236L)

  expect_identical(dim(neighbors), c(4L, 3L))
  expect_type(neighbors, "integer")
  expect_identical(neighbors[, 1L], 0:3)
  expect_true(all(apply(neighbors, 1L, anyDuplicated) == 0L))
  for (cell_index in seq_len(nrow(embedding))) {
    offsets <- sweep(
      embedding[neighbors[cell_index, ] + 1L, , drop = FALSE],
      2L, embedding[cell_index, ], "-"
    )
    expect_equal(rowSums(offsets^2), c(0, 1, 1))
  }
})

test_that("fixed KNN clamps k_max to n_cells - 1", {
  embedding <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 2),
    c(3, 0),
    c(0, 5)
  )

  neighbors <- bcmp:::get_knn_neighbor_once(embedding, k_max = 99L, seed = 236L)
  expect_identical(dim(neighbors), c(5L, 4L))
  expect_true(all(neighbors >= 0L & neighbors < 5L))
})

test_that("fixed KNN validates its direct core inputs", {
  expect_error(
    bcmp:::get_knn_neighbor_once(matrix(1:4, nrow = 2), k_max = 1L),
    "At least 3 cells"
  )
  expect_error(
    bcmp:::get_knn_neighbor_once(matrix(numeric(), nrow = 3, ncol = 0), k_max = 1L),
    "at least one dimension"
  )
  expect_error(
    bcmp:::get_knn_neighbor_once(matrix(1:6, nrow = 3), k_max = 0L),
    "at least 1"
  )
})


test_that("KNN rejects magnitudes unsafe for float32 without changing safe inputs", {
  embedding <- rbind(c(0, 0), c(1, 0), c(0, 2), c(3, 0), c(0, 5))
  expected <- bcmp:::get_knn_neighbor_once(embedding, k_max = 3L)
  for (sparse in c(FALSE, TRUE)) {
    input <- if (sparse) Matrix::Matrix(embedding, sparse = TRUE) else embedding
    expect_identical(
      bcmp:::get_knn_neighbor_once(input * 1e10, k_max = 3L), expected
    )
    # 1e20 fits float32 but its square does not; 1e40 overflows conversion.
    for (scale in c(1e20, 1e40)) {
      expect_error(
        bcmp_embedding(input * scale, c("a", "b", "a", "b", "a"), k_max = 3L,
                       verbose = FALSE),
        "magnitude.*float32"
      )
    }
  }
})

test_that("sparse embeddings reach Annoy without global densification", {
  embedding <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 2),
    c(3, 0),
    c(0, 5)
  )
  sparse_embedding <- Matrix::Matrix(embedding, sparse = TRUE)

  expect_identical(
    bcmp:::get_knn_neighbor_once(sparse_embedding, k_max = 3L, seed = 236L),
    bcmp:::get_knn_neighbor_once(embedding, k_max = 3L, seed = 236L)
  )
})
