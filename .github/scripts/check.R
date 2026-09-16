source(".github/scripts/common.R")
for (package in c("testthat", "jsonlite", "SeuratObject", "SingleCellExperiment",
                  "SummarizedExperiment", "uwot")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing full-check dependency: ", package)
}
Sys.setenv(BCMP_REQUIRE_FULL_REGRESSION = "true", "_R_CHECK_FORCE_SUGGESTS_" = "true")
record_environment()
tarball <- build_tarball()
result <- rcmdcheck::rcmdcheck(
  tarball, args = "--as-cran", check_dir = file.path(output, "check"),
  error_on = "warning", timeout = 2400
)
writeLines(result$notes, file.path(output, "notes.txt"))
# Only the explicitly reviewed prerelease metadata NOTE is accepted.
# Unexpected platform, documentation, dependency or compiled-code NOTEs fail.
for (note in result$notes) {
  lines <- trimws(strsplit(note, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  known <- grepl("^checking CRAN incoming feasibility [.][.][.] *(\\[[^]]+\\] *)?NOTE$", lines) |
    grepl("^Maintainer: ", lines) |
    lines == "New submission" |
    grepl("^Size of tarball: [0-9]+ bytes$", lines)
  if (!all(known)) stop("Unreviewed R CMD check NOTE; see notes.txt.")
}
logs <- list.files(file.path(output, "check"), pattern = "^testthat\\.Rout$",
                   recursive = TRUE, full.names = TRUE)
if (length(logs) != 1L) stop("Expected a completed installed-package test log.")
summary <- grep("\\[ FAIL ", readLines(logs, warn = FALSE), value = TRUE)
if (length(summary) != 1L ||
    !grepl("FAIL 0 \\| WARN 0 \\| SKIP 0 \\| PASS [1-9][0-9]*", summary)) {
  stop("Full check must complete all tests without failure, warning or skip.")
}
writeLines(summary, file.path(output, "test-summary.txt"))
