source(".github/scripts/common.R")
tarball <- build_tarball()
install_tarball(tarball)
library(testthat)
record_environment()
# Start with correctness diagnostics, not a repository-wide formatting rewrite.
linters <- lintr::linters_with_tags("correctness")
paths <- c(list.files("R", pattern = "\\.R$", full.names = TRUE),
           list.files("tests", pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
           list.files(".github/scripts", pattern = "\\.R$", full.names = TRUE))
lints <- do.call(c, lapply(paths, lintr::lint, linters = linters))
capture.output(print(lints), file = file.path(output, "lintr.txt"))
print(lints)
if (length(lints)) stop("lintr diagnostics must be reviewed.")
