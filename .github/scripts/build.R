source(".github/scripts/common.R")
record_environment()
tarball <- build_tarball()
install_tarball(tarball)
setwd(output)
library(bcmp)
stopifnot(is.function(bcmp::bcmp), is.function(bcmp::bcmp_embedding))
# S4 registration adds namespace metadata, not another public function.
x <- matrix(seq_len(120), nrow = 40)
result <- bcmp_embedding(x, rep(c("a", "b"), each = 20),
                         k_max = 3L, verbose = FALSE)
stopifnot(inherits(result, "EmbeddingResult"), length(result$labels) == 40L)
stopifnot(file.exists(system.file("extdata", "diabetic_kidney_lite_seurat.rds",
                                  package = "bcmp", mustWork = TRUE)))
record_environment()
