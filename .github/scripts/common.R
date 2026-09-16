# Shared helpers for isolated CI runs; no package code is sourced from checkout.
# Limit test workers to two and keep BLAS/OpenMP kernels serial.
Sys.setenv(
  RCPP_PARALLEL_NUM_THREADS = "2",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output <- Sys.getenv("BCMP_CI_OUTPUT")
if (!nzchar(output)) stop("Set BCMP_CI_OUTPUT to an independent output directory.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
if (startsWith(paste0(output, "/"), paste0(root, "/"))) {
  stop("BCMP_CI_OUTPUT must be outside the package checkout.")
}

record_environment <- function() {
  capture.output({
    cat("Commit:", Sys.getenv("GITHUB_SHA", unset = "local validation"), "\n")
    print(.libPaths())
    print(installed.packages()[, c("Package", "Version", "LibPath"), drop = FALSE])
    print(sessionInfo())
  }, file = file.path(output, "environment.txt"))
}

run_r <- function(args, log) {
  status <- system2(file.path(R.home("bin"), "R"),
                    c("CMD", vapply(args, shQuote, character(1))),
                    stdout = file.path(output, log), stderr = file.path(output, log))
  if (status != 0L) {
    cat(readLines(file.path(output, log), warn = FALSE), sep = "\n")
    stop("R CMD failed; see ", log)
  }
}

build_tarball <- function() {
  previous <- setwd(output)
  on.exit(setwd(previous))
  run_r(c("build", root), "build.log")
  version <- read.dcf(file.path(root, "DESCRIPTION"))[1L, "Version"]
  tarball <- file.path(output, paste0("bcmp_", version, ".tar.gz"))
  stopifnot(file.exists(tarball))
  writeLines(paste(digest::digest(file = tarball, algo = "sha256"), basename(tarball)),
             file.path(output, "SHA256SUMS"))
  tarball
}

install_tarball <- function(tarball) {
  lib <- file.path(output, "library")
  dir.create(lib)
  run_r(c("INSTALL", "--install-tests", paste0("--library=", lib), tarball),
        "install.log")
  .libPaths(c(lib, .libPaths()))
  stopifnot(normalizePath(find.package("bcmp")) ==
              normalizePath(file.path(lib, "bcmp")))
  lib
}

check_profile <- function(profile, umap) {
  allowed <- switch(profile, core = character(), seurat = "SeuratObject",
                    sce = "SingleCellExperiment",
                    both = c("SeuratObject", "SingleCellExperiment"),
                    stop("Unknown dependency profile."))
  if (umap) allowed <- c(allowed, "uwot")
  for (package in c("SeuratObject", "SingleCellExperiment", "uwot")) {
    installed <- length(find.package(package, quiet = TRUE)) > 0L
    if (package %in% allowed) {
      if (!installed || !requireNamespace(package, quietly = TRUE)) {
        stop("Missing required profile dependency: ", package)
      }
    } else if (installed) {
      stop("Dependency isolation violated: ", package, " is installed.")
    }
  }
  record_environment()
}
