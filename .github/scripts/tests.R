source(".github/scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[2L] %in% c("with-umap", "without-umap"))
profile <- args[1L]
umap <- identical(args[2L], "with-umap")
check_profile(profile, umap)
Sys.setenv(BCMP_REQUIRE_FULL_REGRESSION = "false")
tarball <- build_tarball()
install_tarball(tarball)
setwd(output)
library(bcmp)
# Full golden regression runs only in R CMD check and explicit coverage runs.
filter <- "^(cluster|knn|mixing|preprocess|search|snn|validation|method-registration|object-adapters|object-registration-lifecycle)$"
results <- testthat::test_dir(
  system.file("tests", "testthat", package = "bcmp", mustWork = TRUE),
  filter = filter, package = "bcmp", load_package = "installed",
  reporter = "summary", stop_on_failure = TRUE, stop_on_warning = TRUE
)
summary <- as.data.frame(results)
write.csv(summary[, c("file", "test", "failed", "skipped", "error", "warning", "passed")],
          file.path(output, "tests.csv"), row.names = FALSE)
if (profile == "both" && umap && any(summary$skipped)) {
  stop("Full optional-dependency profile must not skip tests.")
}
record_environment()
