test_that("log1p normalization matches the frozen Python row-scaling rule", {
  counts <- rbind(
    c(1, 0, 3),
    c(0, 2, 0),
    c(0, 0, 0)
  )
  observed <- bcmp:::normalize_log1p_counts(counts)
  expected <- rbind(
    c(log1p(2500), 0, log1p(7500)),
    c(0, log1p(10000), 0),
    c(0, 0, 0)
  )

  expect_equal(observed, expected, tolerance = 1e-14)
})

test_that("log1p normalization preserves sparse structure until needed", {
  counts <- Matrix::Matrix(
    rbind(
      c(1, 0, 3),
      c(0, 2, 0),
      c(0, 0, 0)
    ),
    sparse = TRUE
  )
  observed <- bcmp:::normalize_log1p_counts(counts)

  expect_s4_class(observed, "dgCMatrix")
  expect_equal(as.matrix(observed), bcmp:::normalize_log1p_counts(as.matrix(counts)))
  expect_equal(length(observed@x), 3L)
})

test_that("raw count validation rejects non-finite and negative values", {
  expect_error(
    bcmp:::normalize_log1p_counts(matrix(c(1, -1), nrow = 1L)),
    "non-negative"
  )
  expect_error(
    bcmp:::normalize_log1p_counts(matrix(c(1, Inf), nrow = 1L)),
    "finite"
  )
  expect_error(
    bcmp:::normalize_log1p_counts(matrix(1, nrow = 1L), target_sum = 0),
    "positive"
  )
})

test_that("Seurat VST batch consensus matches the frozen Python fixture", {
  n_cells <- 30L
  n_genes <- 50L
  cell_index <- 0:(n_cells - 1L)
  gene_index <- 0:(n_genes - 1L)
  counts <- outer(cell_index, gene_index, function(i, j) (i * 13 + j * 7 + i * j * 3) %% 17)
  counts[outer(cell_index, gene_index, function(i, j) (i + 2 * j) %% 5 == 0)] <- 0
  batch_labels <- c(rep("a", 15L), rep("b", 15L))
  var_names <- paste0("g", gene_index)

  observed <- bcmp:::select_highly_variable_genes(
    Matrix::Matrix(counts, sparse = TRUE), batch_labels, var_names, 12L
  )
  expected <- c("g39", "g32", "g8", "g15", "g31", "g18", "g2", "g14", "g9", "g17", "g33", "g40")

  expect_identical(observed, expected)
})

test_that("Seurat VST consensus uses frequency, median rank, then feature name", {
  observed <- bcmp:::seurat_v5_consensus_features(
    list(
      counts.a = c("g2", "g1", "g0"),
      counts.b = c("g1", "g2", "g3")
    ),
    common_features = c("g0", "g1", "g2", "g3")
  )

  expect_identical(observed, c("g1", "g2", "g0", "g3"))
})

test_that("duplicate LOESS predictors separate ties and preserve unique values", {
  predictor <- c(0, 0, 1, 2, 2, 3)
  observed <- bcmp:::stabilize_duplicate_loess_predictor(predictor)

  expect_identical(observed[c(3L, 6L)], predictor[c(3L, 6L)])
  expect_true(all(diff(sort(observed)) > 0))
  expect_identical(
    bcmp:::stabilize_duplicate_loess_predictor(predictor), observed
  )
})

test_that("a single LOESS predictor is unchanged", {
  expect_identical(bcmp:::stabilize_duplicate_loess_predictor(1), 1)
})

test_that("Seurat VST consensus handles empty collections and feature overlap", {
  expect_identical(
    bcmp:::seurat_v5_consensus_features(list(), "g1"), character()
  )
  expect_identical(
    bcmp:::seurat_v5_consensus_features(list(counts.a = "g1"), character()),
    character()
  )
  expect_identical(
    bcmp:::seurat_v5_consensus_features(list(counts.a = "g1"), "g2"),
    character()
  )
})

test_that("Seurat VST rejects a batch with fewer than three cells", {
  expect_error(
    bcmp:::select_highly_variable_genes(
      matrix(seq_len(15), nrow = 5L),
      c("tiny", "tiny", "a", "a", "a"),
      c("g0", "g1", "g2"),
      2L
    ),
    "at least 3 cells.*tiny=2"
  )
})

test_that("Seurat-style scaling centers, sample-scales, and upper-clips dense values", {
  matrix <- rbind(
    c(0, 1),
    c(1, 1),
    c(2, 1)
  )
  observed <- bcmp:::seurat_style_scale_dense(matrix, scale_max = 0.5)
  expected <- cbind(c(-1, 0, 0.5), c(0, 0, 0))

  expect_equal(observed, expected, tolerance = 1e-14)
})

test_that("full-SVD PCA preserves the Python-equivalent PCA geometry", {
  scaled <- rbind(
    c(-1, -1),
    c(-1, 1),
    c(1, -1),
    c(1, 1)
  )
  observed <- bcmp:::run_bcmp_pca(scaled, 2L)

  expect_equal(dim(observed), c(4L, 2L))
  expect_true(all(is.finite(observed)))
  expect_equal(as.matrix(stats::dist(observed)), as.matrix(stats::dist(scaled)), tolerance = 1e-12)
  expect_error(bcmp:::run_bcmp_pca(scaled, 3L), "too few PCs")
})

test_that("standard raw-count preprocessing composes the verified Python-aligned stages", {
  n_cells <- 30L
  n_genes <- 50L
  cell_index <- 0:(n_cells - 1L)
  gene_index <- 0:(n_genes - 1L)
  counts <- outer(cell_index, gene_index, function(i, j) (i * 13 + j * 7 + i * j * 3) %% 17)
  counts[outer(cell_index, gene_index, function(i, j) (i + 2 * j) %% 5 == 0)] <- 0

  result <- bcmp:::run_standard_bcmp_preprocess_from_counts(
    raw_counts = Matrix::Matrix(counts, sparse = TRUE),
    batch_labels = c(rep("a", 15L), rep("b", 15L)),
    var_names = paste0("g", gene_index),
    partition_n_hvg = 12L,
    partition_n_pcs = 5L
  )

  expect_identical(
    result$hvg_genes,
    c("g39", "g32", "g8", "g15", "g31", "g18", "g2", "g14", "g9", "g17", "g33", "g40")
  )
  expect_equal(dim(result$pca_embeddings), c(30L, 5L))
  expect_null(result$umap_embeddings)
})

test_that("standard raw-count preprocessing validates the requested PCA dimensions", {
  cell_index <- 0:29
  gene_index <- 0:49
  counts <- outer(cell_index, gene_index, function(i, j) (i * 13 + j * 7 + i * j * 3) %% 17)
  counts[outer(cell_index, gene_index, function(i, j) (i + 2 * j) %% 5 == 0)] <- 0
  expect_error(
    bcmp:::run_standard_bcmp_preprocess_from_counts(
      counts,
      batch_labels = c(rep("a", 15L), rep("b", 15L)),
      var_names = paste0("g", gene_index),
      partition_n_hvg = 2L,
      partition_n_pcs = 3L
    ),
    "only 2 HVGs"
  )
})

test_that("UMAP rejects invalid matrices and too few cells before dependency checks", {
  expect_error(bcmp:::run_bcmp_umap2d(matrix(numeric(), 3L, 0L)),
               "at least 3 cells and 1 PCA dimension")
  expect_error(bcmp:::run_bcmp_umap2d(matrix(NA_real_, 3L, 1L)), "finite")
  expect_error(bcmp:::run_bcmp_umap2d(matrix("invalid", 3L, 1L)), "numeric")
  expect_error(bcmp:::run_bcmp_umap2d(matrix(1, 2L, 1L)), "at least 3 cells")
})

test_that("Seurat-style debug UMAP is seeded, two-dimensional, and RNG-local", {
  bcmp_test_require("uwot")
  index <- seq_len(24L)
  pca_embeddings <- cbind(
    sin(index / 3),
    cos(index / 5),
    index / 10,
    sin(index / 7),
    cos(index / 11)
  )

  set.seed(739L)
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  first <- bcmp:::run_bcmp_umap2d(pca_embeddings, seed = 236L)
  after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  second <- bcmp:::run_bcmp_umap2d(pca_embeddings, seed = 236L)

  expect_identical(after, before)
  expect_equal(first, second, tolerance = 1e-12)
  expect_equal(dim(first), c(24L, 2L))
  expect_identical(colnames(first), c("UMAP_1", "UMAP_2"))
})

test_that("debug preprocessing opt-in adds a Seurat-style UMAP", {
  bcmp_test_require("uwot")
  n_cells <- 30L
  n_genes <- 50L
  cell_index <- 0:(n_cells - 1L)
  gene_index <- 0:(n_genes - 1L)
  counts <- outer(cell_index, gene_index, function(i, j) (i * 13 + j * 7 + i * j * 3) %% 17)
  counts[outer(cell_index, gene_index, function(i, j) (i + 2 * j) %% 5 == 0)] <- 0

  result <- bcmp:::run_standard_bcmp_preprocess_from_counts(
    raw_counts = Matrix::Matrix(counts, sparse = TRUE),
    batch_labels = c(rep("a", 15L), rep("b", 15L)),
    var_names = paste0("g", gene_index),
    partition_n_hvg = 12L,
    partition_n_pcs = 5L,
    compute_umap = TRUE,
    seed = 236L
  )

  expect_equal(dim(result$umap_embeddings), c(30L, 2L))
  expect_identical(colnames(result$umap_embeddings), c("UMAP_1", "UMAP_2"))
})

test_that("raw-count composition reuses the embedding result contract", {
  n_cells <- 30L
  n_genes <- 50L
  cell_index <- 0:(n_cells - 1L)
  gene_index <- 0:(n_genes - 1L)
  counts <- outer(cell_index, gene_index, function(i, j) (i * 13 + j * 7 + i * j * 3) %% 17)
  counts[outer(cell_index, gene_index, function(i, j) (i + 2 * j) %% 5 == 0)] <- 0
  batches <- c(rep("a", 15L), rep("b", 15L))
  cell_ids <- paste0("cell_", cell_index)
  var_names <- paste0("g", gene_index)

  composed <- bcmp:::run_bcmp_from_counts(
    raw_counts = Matrix::Matrix(counts, sparse = TRUE),
    batch_labels = batches,
    var_names = var_names,
    cell_ids = cell_ids,
    partition_n_hvg = 12L,
    partition_n_pcs = 5L,
    k_min = 3L,
    k_max = 4L,
    min_batch_coverage = 2L,
    min_cells_in_domain_per_batch = 1L,
    seed = 236L,
    verbose = FALSE
  )
  direct <- bcmp::bcmp_embedding(
    embedding = composed$preprocess$pca_embeddings,
    batch_labels = batches,
    cell_ids = cell_ids,
    k_min = 3L,
    k_max = 4L,
    min_batch_coverage = 2L,
    min_cells_in_domain_per_batch = 1L,
    seed = 236L,
    verbose = FALSE
  )

  expect_s3_class(composed$core, "EmbeddingResult")
  expect_identical(composed$core$labels, direct$labels)
  expect_equal(composed$core$selection, direct$selection)
  expect_identical(composed$core$search_trace, direct$search_trace)
  expect_null(composed$preprocess$umap_embeddings)
})
