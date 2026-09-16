test_that("public entry points retain their default parameter contract", {
  # Deliberately independent of parameters.rds, which explicit regression calls use.
  common <- list(
    k_min = 3L, k_max = NULL, batch_frac_threshold_for_k_max = 0.01,
    min_batch_coverage = NULL, min_cells_in_domain_per_batch = 5L,
    max_underrepresentation_fold = 10L, max_residual_cell_frac = 0.05,
    seed = 236L, output_level = "standard", verbose = TRUE
  )
  embedding <- c(list(cell_ids = NULL), common)
  object <- c(list(batch_key = "batch", assay = NULL, layer = NULL,
                   partition_n_pcs = 30L, partition_n_hvg = 2000L),
              common[1:7], list(exclude_cells = NULL), common[8:10])
  expect_identical(as.list(formals(bcmp_embedding))[-c(1L, 2L)], embedding)
  expect_identical(as.list(formals(bcmp))[-1L], object)
})

test_that("embedding-only loading does not load object or UMAP dependencies", {
  run_bcmp_registration_process(c(
    'loadNamespace("bcmp")',
    'stopifnot(!any(c("SeuratObject", "SingleCellExperiment", "uwot") %in% loadedNamespaces()))',
    'x <- matrix(seq_len(120), nrow = 40)',
    'b <- rep(c("a", "b"), each = 20)',
    'result <- bcmp::bcmp_embedding(x, b, k_max = 3L, output_level = "debug", verbose = FALSE)',
    'stopifnot(inherits(result, "EmbeddingDebugResult"))',
    'stopifnot(!any(c("SeuratObject", "SingleCellExperiment", "uwot") %in% loadedNamespaces()))'
  ))
})

test_that("unloading and reloading bcmp cleans up its hooks and restores methods", {
  bcmp_test_require("SeuratObject")
  run_bcmp_registration_process(c(
    'loadNamespace("SeuratObject")',
    'event <- packageEvent("SeuratObject", "onLoad")',
    'sentinel <- function(...) invisible(NULL)',
    'setHook(event, sentinel, action = "append")',
    'baseline <- getHook(event)',
    'loadNamespace("bcmp")',
    'stopifnot(length(getHook(event)) == length(baseline) + 1L)',
    'unloadNamespace("bcmp")',
    'stopifnot(identical(getHook(event), baseline))',
    'loadNamespace("bcmp")',
    'stopifnot(length(getHook(event)) == length(baseline) + 1L)',
    'stopifnot(!is.null(methods::selectMethod(bcmp::bcmp, "Seurat", useInherited = FALSE, optional = TRUE)))',
    'unloadNamespace("bcmp")',
    'stopifnot(identical(getHook(event), baseline))'
  ))
})
