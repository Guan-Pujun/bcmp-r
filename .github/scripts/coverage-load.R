# CI-only adaptation of covr's temporary installation, never package sources.
check_load_coverage_support <- function() {
  version <- as.character(utils::packageVersion("covr"))
  if (!identical(version, "3.6.5")) {
    stop("Early-load coverage is validated with covr 3.6.5; review before using ", version)
  }
}

with_early_load_coverage <- function(install_path, code) {
  check_load_coverage_support()
  expected_library <- normalizePath(install_path, winslash = "/", mustWork = FALSE)
  namespace <- asNamespace("covr")
  if (inherits(get("add_hooks", namespace), "functionWithTrace")) {
    stop("covr::add_hooks is already traced; refusing to replace another tracer.")
  }
  rewrite_loader <- function(pkg_name, lib) {
    if (!identical(pkg_name, "bcmp") ||
        !identical(normalizePath(lib, winslash = "/", mustWork = TRUE), expected_library)) {
      stop("Early-load coverage may only modify the designated temporary bcmp installation.")
    }
    loader <- file.path(lib, pkg_name, "R", pkg_name)
    lines <- readLines(loader, warn = FALSE)
    late <- 'setHook(packageEvent(pkg, "onLoad"), function(...) covr:::trace_environment(ns))'
    position <- which(lines == late)
    lazy_load <- grep("lazyLoad(dbbase, ns,", lines, fixed = TRUE)
    if (length(position) != 1L || length(lazy_load) != 1L || position <= lazy_load) {
      stop("Unexpected covr loader layout; early-load coverage was not applied.")
    }
    # R still invokes .onLoad itself. Instrument after lazyLoad but before it runs.
    lines[position] <- "covr:::trace_environment(ns)"
    writeLines(lines, loader)
  }
  trace("add_hooks", where = namespace, print = FALSE,
        exit = substitute(rewrite(pkg_name, lib), list(rewrite = rewrite_loader)))
  on.exit(untrace("add_hooks", where = namespace), add = TRUE)
  force(code)
}

verify_load_coverage <- function(coverage) {
  tally <- covr::tally_coverage(coverage)
  loading <- tally[tally$filename == "R/zzz.R" & tally$functions == ".onLoad", ]
  if (!nrow(loading) || anyNA(loading$value) || any(loading$value <= 0)) {
    stop("Real .onLoad execution was not fully captured; refusing an incomplete coverage report.")
  }
}
