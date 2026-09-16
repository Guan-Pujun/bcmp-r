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
# Private review cannot resolve the formal repository and its unpublished Pages.
# Derive this exception from the GitHub event, never from package metadata.
private_review <- FALSE
event_path <- Sys.getenv("GITHUB_EVENT_PATH")
if (identical(Sys.getenv("GITHUB_ACTIONS"), "true") && file.exists(event_path)) {
  repository <- jsonlite::fromJSON(event_path)$repository
  private_review <- identical(repository$private, TRUE) &&
    identical(repository$full_name, "Guan-Pujun/bcmp-r") &&
    identical(Sys.getenv("GITHUB_REPOSITORY"), repository$full_name)
}
# All other platform, documentation, dependency or compiled-code NOTEs fail.
for (note in result$notes) {
  lines <- trimws(strsplit(note, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  known <- grepl("^checking CRAN incoming feasibility [.][.][.] *(\\[[^]]+\\] *)?NOTE$", lines) |
    grepl("^Maintainer: ", lines) |
    lines == "New submission" |
    grepl("^Size of tarball: [0-9]+ bytes$", lines)
  if (private_review && grepl("^checking CRAN incoming feasibility ", lines[1L])) {
    expected_urls <- c(
      "https://github.com/Guan-Pujun/bcmp-r",
      "https://github.com/Guan-Pujun/bcmp-r/issues",
      "https://guan-pujun.github.io/bcmp-r/"
    )
    accepted_url <- FALSE
    for (url in expected_urls) {
      block <- c(paste("URL:", url), "From: DESCRIPTION",
                 "Status: 404", "Message: Not Found")
      for (start in which(lines == block[1L])) {
        positions <- start + seq_along(block) - 1L
        if (identical(lines[positions], block)) {
          known[positions] <- TRUE
          accepted_url <- TRUE
          message("Private review only: accepting DESCRIPTION URL 404: ", url)
        }
      }
    }
    if (accepted_url) {
      known[lines == "Found the following (possibly) invalid URLs:"] <- TRUE
    }
  }
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
