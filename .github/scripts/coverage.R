source(".github/scripts/common.R")
for (package in c("covr", "xml2", "digest", "testthat", "jsonlite",
                  "SeuratObject", "SingleCellExperiment", "SummarizedExperiment", "uwot")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing coverage dependency: ", package)
}
source(".github/scripts/coverage-load.R")
check_load_coverage_support()
# Exercise the real exporter before the full run, using one trivial R expression.
probe_env <- new.env()
eval(parse(text = "coverage_probe <- function(x) {\n x + 1\n}",
           keep.source = TRUE), envir = probe_env)
probe <- covr::function_coverage(
  "coverage_probe", quote(coverage_probe(1)), env = probe_env, enc = probe_env
)
probe_file <- file.path(output, "coverage-preflight.xml")
covr::to_cobertura(probe, filename = probe_file)
probe_xml <- xml2::read_xml(probe_file)
stopifnot(identical(xml2::xml_name(probe_xml), "coverage"),
          identical(xml2::xml_attr(probe_xml, "lines-valid"), "1"),
          identical(xml2::xml_attr(probe_xml, "lines-covered"), "1"))
if ("--preflight-only" %in% commandArgs(trailingOnly = TRUE)) {
  quit(save = "no", status = 0L)
}
Sys.setenv(BCMP_REQUIRE_FULL_REGRESSION = "true")
record_environment()
tarball <- build_tarball()
source_dir <- file.path(output, "coverage-source")
dir.create(source_dir)
utils::untar(tarball, exdir = source_dir)
coverage <- with_early_load_coverage(file.path(output, "coverage-library"),
  covr::package_coverage(
    path = file.path(source_dir, "bcmp"), type = "tests",
    install_path = file.path(output, "coverage-library"),
    quiet = FALSE, clean = FALSE
  )
)
# Preserve the result before reporting, so export failures do not lose the run.
saveRDS(coverage, file.path(output, "coverage.rds"))
verify_load_coverage(coverage)
capture.output({
  cat("Overall coverage:", covr::percent_coverage(coverage), "\n")
  print(covr::tally_coverage(coverage))
}, file = file.path(output, "coverage-summary.txt"))
# Include R and compiled-code coverage; do not exclude important modules.
xml_file <- file.path(output, "coverage.xml")
covr::to_cobertura(coverage, filename = xml_file)
xml <- xml2::read_xml(xml_file)
stopifnot(identical(xml2::xml_name(xml), "coverage"),
          as.integer(xml2::xml_attr(xml, "lines-valid")) > 0L,
          length(xml2::xml_find_all(xml, ".//class")) > 0L)
# Match the Python CI coverage gate without changing the regression rules.
percent <- covr::percent_coverage(coverage)
if (!is.finite(percent) || percent < 90) {
  stop("Coverage must be at least 90%; see coverage-summary.txt.")
}
