# Keep optional S4 registrations outside the sealed package namespace.
bcmp_object_hooks <- new.env(parent = emptyenv())

register_bcmp_object_method <- function(package) {
  if (!isNamespaceLoaded(package)) return(invisible(FALSE))
  adapter <- switch(
    package,
    SeuratObject = list(class = "Seurat", definition = bcmp_seurat),
    SingleCellExperiment = list(class = "SingleCellExperiment", definition = bcmp_sce)
  )
  if (is.null(adapter)) stop("Unknown BCMP object adapter: ", package, call. = FALSE)
  methods::setMethod(
    "bcmp", c(object = adapter$class), definition = adapter$definition,
    where = bcmp_object_hooks$registry
  )
  bcmp_object_hooks$registered[[package]] <- adapter$class
  invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
  bcmp_object_hooks$registry <- new.env(parent = asNamespace(pkgname))
  bcmp_object_hooks$registered <- list()
  bcmp_object_hooks$callbacks <- list()
  for (package in c("SeuratObject", "SingleCellExperiment")) {
    callback <- local({
      adapter_package <- package
      function(...) register_bcmp_object_method(adapter_package)
    })
    bcmp_object_hooks$callbacks[[package]] <- callback
    setHook(packageEvent(package, "onLoad"), callback, action = "append")
    callback()
  }
}

.onUnload <- function(libpath) {
  for (package in names(bcmp_object_hooks$callbacks)) {
    event <- packageEvent(package, "onLoad")
    callback <- bcmp_object_hooks$callbacks[[package]]
    remaining <- Filter(function(hook) !identical(hook, callback), getHook(event))
    setHook(event, remaining, action = "replace")
  }
  for (class in bcmp_object_hooks$registered) {
    methods::removeMethod("bcmp", class, where = bcmp_object_hooks$registry)
  }
  bcmp_object_hooks$callbacks <- list()
  bcmp_object_hooks$registered <- list()
}
