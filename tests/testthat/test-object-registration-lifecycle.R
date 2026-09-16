for (adapter_package in c("SeuratObject", "SingleCellExperiment")) {
  local({
    package <- adapter_package
    class <- if (package == "SeuratObject") "Seurat" else "SingleCellExperiment"
    other <- if (package == "SeuratObject") "SingleCellExperiment" else "SeuratObject"
    construct <- if (package == "SeuratObject") {
      'object <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts)); object$batch <- batches'
    } else {
      'object <- SingleCellExperiment::SingleCellExperiment(list(counts = counts)); SummarizedExperiment::colData(object)$batch <- batches'
    }
    setup <- c(
      'set.seed(1107)',
      'counts <- Matrix::Matrix(matrix(rpois(60 * 104, 5), 60, dimnames = list(paste0("g", 1:60), paste0("c", 1:104))), sparse = TRUE)',
      'batches <- rep(c("a", "b"), each = 52)',
      construct
    )
    check <- c(
      'result <- bcmp::bcmp(object, partition_n_pcs = 10L, partition_n_hvg = 20L, k_max = 3L, verbose = FALSE)',
      'stopifnot(inherits(result, "ObjectResult"), result$selection$selected_k == 3L)',
      'stopifnot(!isNamespaceLoaded("uwot"))'
    )

    test_that(paste(package, "runs the full workflow in either loading order"), {
      bcmp_test_require(package)
      for (dependency_first in c(FALSE, TRUE)) {
        order <- if (dependency_first) c(package, "bcmp") else c("bcmp", package)
        run_bcmp_registration_process(c(
          sprintf('loadNamespace("%s")', order),
          setup, check,
          sprintf('stopifnot(!isNamespaceLoaded("%s"))', other),
          sprintf('before <- methods::selectMethod(bcmp::bcmp, "%s", useInherited = FALSE)', class),
          'stopifnot(identical(formals(before), formals(bcmp::bcmp)))',
          sprintf('bcmp:::register_bcmp_object_method("%s")', package),
          sprintf('after <- methods::selectMethod(bcmp::bcmp, "%s", useInherited = FALSE)', class),
          # covr replaces method bodies with counters after namespace loading.
          # Preserve identity checks normally and verify real output in both modes.
          'if (!identical(Sys.getenv("R_COVR"), "true")) stopifnot(identical(before, after))',
          'stopifnot(identical(formals(before), formals(after)))',
          'previous_result <- result',
          check,
          'stopifnot(identical(previous_result, result))'
        ))
      }
    })

    test_that(paste(package, "dispatches after dependency reload and for subclasses"), {
      bcmp_test_require(package)
      run_bcmp_registration_process(c(
        'loadNamespace("bcmp")',
        sprintf('loadNamespace("%s")', package),
        sprintf('unloadNamespace("%s")', package),
        sprintf('loadNamespace("%s")', package),
        setup,
        check,
        sprintf('methods::setClass("BcmpRegistrationSubclass", contains = "%s")', class),
        'slots <- setNames(lapply(methods::slotNames(object), function(n) methods::slot(object, n)), methods::slotNames(object))',
        'object <- do.call(methods::new, c(list(Class = "BcmpRegistrationSubclass"), slots))',
        check,
        'unloadNamespace("bcmp")',
        'loadNamespace("bcmp")',
        check
      ))
    })

    test_that(paste(package, "dispatches a serialized object in a fresh session"), {
      bcmp_test_require(package)
      object_file <- tempfile("bcmp-registration-object-", fileext = ".rds")
      path <- paste(deparse(object_file), collapse = "")
      run_bcmp_registration_process(c(setup, sprintf('saveRDS(object, %s)', path)))
      run_bcmp_registration_process(c(
        'loadNamespace("bcmp")',
        sprintf('stopifnot(!isNamespaceLoaded("%s"))', package),
        sprintf('object <- readRDS(%s)', path),
        check
      ))
    })
  })
}
