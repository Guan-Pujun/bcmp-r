# R-only frozen outputs. Updating the fixture requires explicit review.
r_golden_path <- function(file) {
  path <- test_path("fixtures", "r-v1", "diabetic-kidney-lite", file)
  if (!file.exists(path)) stop("Missing required R regression fixture: ", file)
  path
}

# Match Python's cell-aligned, literal-label mismatch count.
r_label_mismatches <- function(observed, expected) {
  stopifnot(length(observed) > 0L, length(observed) == length(expected),
            !anyNA(observed), !anyNA(expected))
  if (!identical(names(observed), names(expected))) {
    stop("Cell IDs and order must match the golden labels.", call. = FALSE)
  }
  as.integer(sum(as.character(observed) != as.character(expected)))
}

expect_r_golden_labels <- function(observed, expected) {
  mismatches <- r_label_mismatches(observed, expected)
  # Integer arithmetic makes exactly 1% fail, including small fixtures.
  expect_lt(mismatches * 100, length(expected))
}

expect_r_domain_summary <- function(observed, expected) {
  # Compare the complete table: do not align away extra or duplicate rows.
  expect_equal(observed, expected, tolerance = 1e-12)
}

expect_r_debug <- function(observed, expected) {
  expect_identical(names(observed), names(expected))
  expect_identical(observed$parameters, expected$parameters)
  expect_identical(observed$search, expected$search)
  if ("workflow" %in% names(expected)) {
    expect_identical(names(observed$workflow), names(expected$workflow))
    expect_setequal(observed$workflow$hvg_genes, expected$workflow$hvg_genes)
    other_fields <- setdiff(names(expected$workflow), "hvg_genes")
    expect_identical(observed$workflow[other_fields],
                     expected$workflow[other_fields])
  }
  expect_identical(names(observed$candidates), names(expected$candidates))
  for (key in names(expected$candidates)) {
    a <- observed$candidates[[key]]
    b <- expected$candidates[[key]]
    expect_identical(names(a), names(b))
    expect_r_golden_labels(a$partition, b$partition)
    expect_r_domain_summary(a$domain_summary, b$domain_summary)
  }
}

expect_r_core <- function(result, expected, labels) {
  expect_identical(class(result), expected$result_class)
  expect_identical(names(result), expected$result_fields)
  expect_equal(result$selection, expected$selection, tolerance = 1e-12)
  expect_equal(result$search_trace, expected$search_trace, tolerance = 1e-12)
  expect_r_golden_labels(labels, expected$labels)
  if ("debug" %in% names(expected)) {
    expect_r_debug(result$debug, expected$debug)
  } else {
    expect_false("debug" %in% names(result))
  }
}

test_that("golden label tolerance matches Python's strict one percent rule", {
  expected <- stats::setNames(rep(c("A", "B", "C"), each = 500L),
                              as.character(seq_len(1500L)))
  observed <- expected
  expect_r_golden_labels(observed, expected)
  observed[seq_len(14L)] <- "B"
  expect_identical(r_label_mismatches(observed, expected), 14L)
  expect_r_golden_labels(observed, expected)
  observed[15L] <- "B"
  expect_identical(r_label_mismatches(observed, expected), 15L)
  expect_failure(expect_r_golden_labels(observed, expected))

  # Python compares label values; it does not optimize a relabeling.
  observed[] <- paste0("renamed_", expected)
  expect_identical(r_label_mismatches(observed, expected), 1500L)
  expect_failure(expect_r_golden_labels(observed, expected))
})

test_that("golden label comparison rejects cell order misalignment", {
  # Identical values ensure only the cell order can trigger this failure.
  expected <- stats::setNames(rep("A", 2L), c("cell_1", "cell_2"))
  expect_error(expect_r_golden_labels(rev(expected), expected),
               "Cell IDs and order must match")
})

test_that("debug workflow compares HVGs as a set", {
  expected <- list(
    parameters = list(), search = list(), candidates = list(),
    workflow = list(hvg_genes = c("gene_a", "gene_b"),
                    excluded_cell_ids = character())
  )
  observed <- expected
  observed$workflow$hvg_genes <- rev(observed$workflow$hvg_genes)
  expect_r_debug(observed, expected)
})

test_that("debug candidate labels use the same tolerance as final labels", {
  expected <- readRDS(r_golden_path("embedding-debug.rds"))$debug
  observed <- expected
  # Swap two cells to preserve domain sizes while changing their assignments.
  for (key in names(observed$candidates)) {
    partition <- observed$candidates[[key]]$partition
    other <- which(partition != partition[[1L]])[1L]
    if (!is.na(other)) {
      partition[c(1L, other)] <- partition[c(other, 1L)]
      observed$candidates[[key]]$partition <- partition
    }
  }
  expect_r_debug(observed, expected)
})

test_that("R regression fixtures and bundled input are present", {
  files <- c(
    "input-manifest.rds", "parameters.rds", "embedding-input.rds",
    "seurat-standard.rds", "seurat-debug.rds",
    "sce-standard.rds", "sce-debug.rds",
    "embedding-standard.rds", "embedding-debug.rds",
    "selection-summary.json", "checksums.csv"
  )
  for (file in files) expect_true(file.exists(r_golden_path(file)))
  input <- system.file("extdata", "diabetic_kidney_lite_seurat.rds",
                       package = "bcmp", mustWork = TRUE)
  expect_identical(unname(tools::md5sum(input)),
                   "769b6bcff84cc8a1ea006b1cc53afde1")
})

for (level in c("standard", "debug")) {
  local({
    output_level <- level
    test_that(paste("R embedding", output_level, "matches its own golden"), {
      input <- readRDS(r_golden_path("embedding-input.rds"))
      args <- readRDS(r_golden_path("parameters.rds"))$embedding
      args$output_level <- output_level
      expected <- readRDS(r_golden_path(
        paste0("embedding-", output_level, ".rds")
      ))
      result <- do.call(bcmp_embedding, c(input, args))
      expect_r_core(result, expected, result$labels)
    })
  })
}

# Object adapters retain raw-count standard regressions; their debug output
# contracts are exercised by the small object-adapter tests.
for (adapter in c("seurat", "sce")) {
  local({
    object_adapter <- adapter
    output_level <- "standard"
    test_that(paste("R", object_adapter, output_level,
                    "matches its own golden"), {
      # The bundled raw-count example is serialized as a Seurat object.
      packages <- "SeuratObject"
      if (object_adapter == "sce") {
        packages <- c(packages, "SingleCellExperiment", "SummarizedExperiment")
      }
      bcmp_test_require(packages)

      manifest <- readRDS(r_golden_path("input-manifest.rds"))
      seurat <- readRDS(system.file(
        "extdata", "diabetic_kidney_lite_seurat.rds",
        package = "bcmp", mustWork = TRUE
      ))
      counts <- SeuratObject::LayerData(seurat, assay = "RNA", layer = "counts")
      cells <- colnames(counts)
      batch <- stats::setNames(as.character(seurat[[]][cells, "batch"]), cells)
      expect_identical(cells, manifest$cells)
      expect_identical(rownames(counts), manifest$genes)
      expect_identical(batch, manifest$batch)
      object <- seurat
      if (object_adapter == "sce") {
        object <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
        SummarizedExperiment::colData(object)$batch <- unname(batch)
      }
      before <- serialize(object, NULL, version = 2)
      args <- readRDS(r_golden_path("parameters.rds"))$object
      args$output_level <- output_level
      expected <- readRDS(r_golden_path(
        paste0(object_adapter, "-", output_level, ".rds")
      ))
      result <- do.call(bcmp, c(list(object = object), args))
      expect_identical(serialize(object, NULL, version = 2), before)
      output <- result$object
      expect_identical(colnames(output), cells)
      if (object_adapter == "seurat") {
        labels <- stats::setNames(as.character(output[[]][cells, "bcmp_domain"]), cells)
        actual_counts <- SeuratObject::LayerData(output, assay = "RNA", layer = "counts")
        actual_batch <- as.character(output[[]][cells, "batch"])
        reductions <- SeuratObject::Reductions(output)
      } else {
        labels <- stats::setNames(as.character(
          SummarizedExperiment::colData(output)$bcmp_domain
        ), cells)
        actual_counts <- SummarizedExperiment::assay(output, "counts")
        actual_batch <- as.character(SummarizedExperiment::colData(output)$batch)
        reductions <- SingleCellExperiment::reducedDimNames(output)
      }
      expect_identical(actual_counts, counts)
      expect_identical(actual_batch, unname(batch))
      expect_r_core(result, expected, labels)
      expect_identical(intersect(c("BCMP_PCA", "BCMP_UMAP"), reductions),
                       expected$bcmp_reduction_names)
    })
  })
}
