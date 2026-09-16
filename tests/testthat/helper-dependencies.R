# Optional runs may skip unavailable ecosystems; designated full runs must fail.
bcmp_test_require <- function(packages) {
  strict <- identical(Sys.getenv("BCMP_REQUIRE_FULL_REGRESSION"), "true")
  for (package in packages) {
    if (strict && !requireNamespace(package, quietly = TRUE)) {
      stop("Full R validation requires package: ", package)
    }
    skip_if_not_installed(package)
  }
}
