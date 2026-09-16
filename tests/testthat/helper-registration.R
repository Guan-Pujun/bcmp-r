run_bcmp_registration_process <- function(code) {
  package_path <- find.package("bcmp")
  installed <- file.exists(file.path(package_path, "Meta", "package.rds"))
  if (!installed && identical(Sys.getenv("BCMP_REQUIRE_FULL_REGRESSION"), "true")) {
    stop("Full R validation requires an installed package for fresh-process tests.")
  }
  skip_if(!installed, "Fresh-process registration tests require an installed package.")
  script <- tempfile("bcmp-registration-", fileext = ".R")
  paths <- unique(c(dirname(package_path), .libPaths()))
  setup <- paste0(".libPaths(", paste(deparse(paths), collapse = "\n"), ")")
  writeLines(c(setup, code), script)
  previous_r_tests <- Sys.getenv("R_TESTS", unset = NA_character_)
  on.exit({
    if (is.na(previous_r_tests)) Sys.unsetenv("R_TESTS")
    else Sys.setenv(R_TESTS = previous_r_tests)
  }, add = TRUE)
  # system2(env=) passes command-line arguments on Windows, not an environment.
  Sys.setenv(R_TESTS = "")
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(script)),
    stdout = TRUE, stderr = TRUE, timeout = 120
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(c(output, attr(output, "errmsg")), collapse = "\n"))
}
