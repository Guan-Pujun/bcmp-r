make_bcmp_adapter_counts <- function(n_cells = 104L, n_genes = 60L) {
  set.seed(1107)
  gene_means <- rep(seq(2, 14, length.out = n_genes), n_cells)
  counts <- matrix(stats::rpois(n_cells * n_genes, lambda = gene_means), nrow = n_genes)
  rownames(counts) <- paste0("gene-", seq_len(n_genes))
  colnames(counts) <- paste0("cell_", seq_len(n_cells))
  Matrix::Matrix(counts, sparse = TRUE)
}

bcmp_adapter_args <- function(...) {
  c(list(partition_n_pcs = 10L, partition_n_hvg = 20L, k_min = 3L,
         k_max = 3L, seed = 236L, verbose = FALSE), list(...))
}

test_that("unsupported bcmp inputs fail explicitly", {
  expect_error(bcmp(list()), "supports only Seurat v5 and SingleCellExperiment")
})

test_that("Seurat v5 standard output preserves input and removes BCMP reductions", {
  bcmp_test_require("SeuratObject")
  counts <- make_bcmp_adapter_counts()
  object <- SeuratObject::CreateSeuratObject(counts, min.cells = 0L, min.features = 0L)
  object$batch <- rep(c("batch_a", "batch_b"), length.out = ncol(object))
  stale <- matrix(0, nrow = ncol(object), ncol = 2L, dimnames = list(colnames(object), c("STALE_1", "STALE_2")))
  object[["BCMP_PCA"]] <- SeuratObject::CreateDimReducObject(stale, assay = "RNA", key = "STALE_")
  result <- do.call(bcmp, c(list(object = object), bcmp_adapter_args()))
  expect_s3_class(result, "ObjectResult")
  expect_false(inherits(result, "ObjectDebugResult"))
  expect_true("bcmp_domain" %in% names(result$object[[]]))
  expect_false("BCMP_PCA" %in% SeuratObject::Reductions(result$object))
  expect_false("BCMP_UMAP" %in% SeuratObject::Reductions(result$object))
  expect_true("BCMP_PCA" %in% SeuratObject::Reductions(object))
  expect_false("bcmp_domain" %in% names(object[[]]))
  expect_error(bcmp(object, layer = "invalid"), "available raw-count layer")
  expect_named(result$selection, c("selected_k", "n_domains", "mixing_status", "n_evaluated_domains", "n_residual_domains", "n_residual_cells", "residual_cell_frac"))
})

test_that("Seurat objects with a legacy Assay receive the supported-version error", {
  bcmp_test_require("SeuratObject")
  assay <- SeuratObject::CreateAssayObject(counts = make_bcmp_adapter_counts())
  object <- SeuratObject::CreateSeuratObject(counts = assay)
  expect_error(bcmp(object, verbose = FALSE), "Seurat v5 Assay5 objects only")
})

test_that("Seurat v5 debug output writes BCMP reductions and excluded labels", {
  bcmp_test_require(c("SeuratObject", "uwot"))
  counts <- make_bcmp_adapter_counts()
  object <- SeuratObject::CreateSeuratObject(counts, min.cells = 0L, min.features = 0L)
  object$batch <- rep(c("batch_a", "batch_b"), length.out = ncol(object))
  excluded <- colnames(object)[c(1L, 3L)]
  result <- do.call(bcmp, c(list(object = object), bcmp_adapter_args(exclude_cells = excluded, output_level = "debug")))
  expect_s3_class(result, "ObjectDebugResult")
  expect_s3_class(result, "ObjectResult")
  expect_true(all(c("BCMP_PCA", "BCMP_UMAP") %in% SeuratObject::Reductions(result$object)))
  expect_true(all(is.na(SeuratObject::Embeddings(result$object, "BCMP_PCA")[excluded, , drop = FALSE])))
  domain_frame <- result$object[["bcmp_domain"]]
  expect_identical(as.character(domain_frame[match(excluded, rownames(domain_frame)), 1L]), rep("Excluded_batch_a", length(excluded)))
  expect_true(all(c("k_min", "k_max", "output_level", "verbose") %in% names(result$debug$parameters)))
})

test_that("SingleCellExperiment standard output preserves input and removes BCMP reductions", {
  bcmp_test_require("SingleCellExperiment")
  counts <- make_bcmp_adapter_counts()
  object <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  SummarizedExperiment::colData(object)$batch <- rep(c("batch_a", "batch_b"), length.out = ncol(object))
  colnames(object) <- colnames(counts)
  SingleCellExperiment::reducedDim(object, "BCMP_PCA") <- matrix(0, nrow = ncol(object), ncol = 2L)
  expect_identical(colnames(object), colnames(counts))
  standard <- do.call(bcmp, c(list(object = object), bcmp_adapter_args()))
  expect_s3_class(standard, "ObjectResult")
  expect_true("bcmp_domain" %in% colnames(SummarizedExperiment::colData(standard$object)))
  expect_false("BCMP_PCA" %in% SingleCellExperiment::reducedDimNames(standard$object))
  expect_true("BCMP_PCA" %in% SingleCellExperiment::reducedDimNames(object))
  expect_error(bcmp(object, layer = "counts"), "only supported for Seurat v5")
})

test_that("SingleCellExperiment debug output writes BCMP reductions and excluded labels", {
  bcmp_test_require(c("SingleCellExperiment", "uwot"))
  counts <- make_bcmp_adapter_counts()
  object <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  SummarizedExperiment::colData(object)$batch <- rep(c("batch_a", "batch_b"), length.out = ncol(object))
  excluded <- colnames(object)[c(1L, 3L)]
  debug <- do.call(bcmp, c(list(object = object), bcmp_adapter_args(exclude_cells = excluded, output_level = "debug")))
  expect_s3_class(debug, "ObjectDebugResult")
  expect_s3_class(debug, "ObjectResult")
  expect_true(all(c("BCMP_PCA", "BCMP_UMAP") %in% SingleCellExperiment::reducedDimNames(debug$object)))
  expect_true(all(is.na(SingleCellExperiment::reducedDim(debug$object, "BCMP_PCA")[excluded, , drop = FALSE])))
  expect_identical(as.character(SummarizedExperiment::colData(debug$object)$bcmp_domain[c(1L, 3L)]), rep("Excluded_batch_a", 2L))
})

for (adapter_package in c("SeuratObject", "SingleCellExperiment")) {
  local({
    package <- adapter_package
    make_object <- function() {
      counts <- make_bcmp_adapter_counts()
      batches <- rep(c("a", "b"), length.out = ncol(counts))
      if (package == "SeuratObject") {
        object <- SeuratObject::CreateSeuratObject(counts)
        object$batch <- batches
      } else {
        object <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
        SummarizedExperiment::colData(object)$batch <- batches
      }
      object
    }

    test_that(paste(package, "rejects an assay not present in the object"), {
      bcmp_test_require(package)
      expect_error(
        bcmp(make_object(), assay = "missing", verbose = FALSE),
        "`assay` must name an assay present"
      )
    })

    test_that(paste(package, "preserves one-PC core results while omitting debug UMAP"), {
      bcmp_test_require(package)
      if (identical(Sys.getenv("BCMP_REQUIRE_FULL_REGRESSION"), "true")) {
        bcmp_test_require("uwot")
      }
      object <- make_object()
      stale <- matrix(0, nrow = ncol(object), ncol = 2L,
                      dimnames = list(colnames(object), c("STALE_1", "STALE_2")))
      if (package == "SeuratObject") {
        object[["BCMP_UMAP"]] <- SeuratObject::CreateDimReducObject(
          stale, assay = "RNA", key = "STALE_")
      } else {
        SingleCellExperiment::reducedDim(object, "BCMP_UMAP") <- stale
      }
      before <- serialize(object, NULL)
      args <- c(list(object = object), bcmp_adapter_args())
      args$partition_n_pcs <- 1L
      standard <- do.call(bcmp, args)
      expect_s3_class(standard, "ObjectResult")
      args$output_level <- "debug"
      if (!requireNamespace("uwot", quietly = TRUE)) {
        expect_error(do.call(bcmp, args), "requires the suggested package `uwot`")
      } else {
        expect_warning(debug <- do.call(bcmp, args),
                       "Skipping BCMP debug UMAP because only one PCA dimension")
        expect_s3_class(debug, "ObjectDebugResult")
        expect_identical(debug$object$bcmp_domain, standard$object$bcmp_domain)
        expect_identical(debug$selection, standard$selection)
        expect_identical(debug$search_trace, standard$search_trace)
        if (package == "SeuratObject") {
          pca <- SeuratObject::Embeddings(debug$object, "BCMP_PCA")
          reductions <- SeuratObject::Reductions(debug$object)
        } else {
          pca <- SingleCellExperiment::reducedDim(debug$object, "BCMP_PCA")
          reductions <- SingleCellExperiment::reducedDimNames(debug$object)
        }
        expect_equal(dim(pca), c(ncol(object), 1L))
        expect_identical(rownames(pca), colnames(object))
        expect_false("BCMP_UMAP" %in% reductions)
      }
      expect_identical(serialize(object, NULL), before)
    })

    test_that(paste(package, "accepts equivalent cell-name and logical exclusions"), {
      bcmp_test_require(package)
      object <- make_object()
      excluded <- colnames(object)[c(1L, 3L)]
      args <- c(list(object = object), bcmp_adapter_args(exclude_cells = excluded))
      named <- do.call(bcmp, args)
      args$exclude_cells <- colnames(object) %in% excluded
      masked <- do.call(bcmp, args)
      expect_length(masked$object$bcmp_domain, ncol(object))
      expect_identical(masked$selection, named$selection)
      expect_identical(masked$search_trace, named$search_trace)
      expect_identical(masked$object$bcmp_domain, named$object$bcmp_domain)
      expect_identical(colnames(masked$object), colnames(object))
    })

    test_that(paste(package, "rejects missing or empty batch metadata"), {
      bcmp_test_require(package)
      object <- make_object()
      expect_error(
        do.call(bcmp, c(list(object = object), bcmp_adapter_args(batch_key = "missing"))),
        "Missing batch metadata column: missing"
      )
      for (label in list(NA_character_, "", NA_real_, NaN)) {
        batches <- if (is.numeric(label)) rep(1, ncol(object)) else object$batch
        batches[1L] <- label
        object$batch <- batches
        expect_error(
          do.call(bcmp, c(list(object = object), bcmp_adapter_args())),
          "Batch metadata column .*batch.*non-missing, non-empty"
        )
      }
    })

    test_that(paste(package, "rejects malformed cell exclusions"), {
      bcmp_test_require(package)
      object <- make_object()
      args <- c(list(object = object), bcmp_adapter_args())
      for (mask in list(TRUE, rep(NA, ncol(object)))) {
        args$exclude_cells <- mask
        expect_error(do.call(bcmp, args), "one non-missing value per cell")
      }
      args$exclude_cells <- 1L
      expect_error(do.call(bcmp, args), "character vector of cell IDs")
      args$exclude_cells <- NA_character_
      expect_error(do.call(bcmp, args), "missing or empty cell IDs")
      args$exclude_cells <- "unknown-cell"
      expect_error(do.call(bcmp, args), "not present in the object: unknown-cell")
    })

    test_that(paste(package, "rejects exclusions leaving fewer than 100 active cells"), {
      bcmp_test_require(package)
      object <- make_object()
      args <- c(
        list(object = object),
        bcmp_adapter_args(exclude_cells = colnames(object)[seq_len(5L)])
      )
      expect_error(
        do.call(bcmp, args),
        "Need at least 100 active cells after applying .*exclude_cells.*got 99"
      )
    })

    test_that(paste(package, "rejects fractional integer parameters before preprocessing"), {
      bcmp_test_require(package)
      object <- make_object()
      for (name in c("partition_n_pcs", "partition_n_hvg", "k_min", "k_max",
                     "min_batch_coverage", "min_cells_in_domain_per_batch",
                     "seed")) {
        args <- list(object = object, verbose = FALSE)
        args[[name]] <- 2.5
        expect_error(do.call(bcmp, args), name)
      }
    })

    test_that(paste(package, "preserves fractional fold settings in debug output"), {
      bcmp_test_require(c(package, "uwot"))
      object <- make_object()
      args <- c(list(object = object), bcmp_adapter_args(
        max_underrepresentation_fold = 3.33, output_level = "debug"))
      result <- do.call(bcmp, args)
      expect_identical(result$debug$parameters$max_underrepresentation_fold, 3.33)
      args$max_underrepresentation_fold <- 0.5
      expect_error(do.call(bcmp, args), "max_underrepresentation_fold")
    })
  })
}
