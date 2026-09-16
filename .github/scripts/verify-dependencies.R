source(".github/scripts/common.R")

# These dependencies re-export the same function. Retained source-reference
# environments can make the copies non-identical across package serialization.
# Install CI dependencies without retained sources; keep sources for BCMP itself.
withCallingHandlers(
  loadNamespace("SummarizedExperiment"),
  warning = function(condition) stop(condition)
)
original <- getExportedValue("S4Arrays", "makeNindexFromArrayViewport")
reexport <- getExportedValue("DelayedArray", "makeNindexFromArrayViewport")
stopifnot(identical(original, reexport))
record_environment()
writeLines("Bioconductor re-exports agree; SummarizedExperiment loads without warnings.",
           file.path(output, "dependency-preflight.txt"))
