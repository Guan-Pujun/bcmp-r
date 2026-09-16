#' Run BCMP from a supported single-cell object
#'
#' Run the complete BCMP workflow from raw counts in a Seurat v5 or
#' SingleCellExperiment object. BCMP selects highly variable genes within
#' batches, normalizes and scales the counts, computes an unintegrated PCA
#' representation, and runs the partition search.
#'
#' @param object (`Seurat` / `SingleCellExperiment`) A Seurat object with an `Assay5`
#'   assay, or a SingleCellExperiment object.
#' @param batch_key (`character`) Name of the batch metadata column in Seurat or
#'   `colData`.
#' @param assay (`character` / `NULL`) Name of the raw-count assay. `NULL` uses the
#'   default Seurat assay or `"counts"` for SingleCellExperiment.
#' @param layer (`character` / `NULL`) Name of the Seurat v5 raw-count layer, e.g.
#'   `"counts"`. `NULL` selects `"counts"`. For SingleCellExperiment, leave this as
#'   `NULL` and select the raw-count assay with `assay`.
#' @param partition_n_pcs (`integer`) Number of principal components used for BCMP.
#' @param partition_n_hvg (`integer`) Number of highly variable genes used for BCMP.
#' @param k_min (`integer`) Lower bound of the representation-scale search.
#' @param k_max (`integer` / `NULL`) Optional upper bound. `NULL` derives it from
#'   batch sizes.
#' @param batch_frac_threshold_for_k_max (`numeric`) Batches at or below this global
#'   cell fraction are ignored during automatic upper-bound derivation.
#' @param min_batch_coverage (`integer` / `NULL`) Minimum number of batches required
#'   in each evaluated domain. `NULL` resolves to `max(2, ceiling(n_batches * 0.5))`.
#' @param min_cells_in_domain_per_batch (`integer`) Minimum effective cells required
#'   for a batch to count as represented in a domain.
#' @param max_underrepresentation_fold (`numeric`) Maximum allowed fold
#'   under-representation relative to global batch proportion.
#' @param max_residual_cell_frac (`numeric`) Maximum fraction of active cells that
#'   may fall in residual domains excluded from the pass/fail decision.
#' @param exclude_cells (`character` / `factor` / `logical` / `NULL`) Cell IDs or a
#'   logical mask in object cell order, where `TRUE` excludes a cell from
#'   partitioning. Excluded cells remain in the returned object as
#'   `Excluded_<batch>`.
#' @param seed (`integer`) Random seed used by preprocessing, partitioning, and debug
#'   UMAP.
#' @param output_level (`character`) `"standard"` (default) or `"debug"`.
#' @param verbose (`logical`) Emit progress records when `TRUE`.
#'
#' @return `"standard"` returns an `ObjectResult` with the copied `object`,
#'   `selection`, and `search_trace`. `"debug"` returns an `ObjectDebugResult`,
#'   which adds `debug`.
#'
#' @aliases bcmp,ANY-method bcmp,Seurat-method bcmp,SingleCellExperiment-method
#' @export
methods::setGeneric(
  "bcmp",
  function(
    object,
    batch_key = "batch",
    assay = NULL,
    layer = NULL,
    partition_n_pcs = 30L,
    partition_n_hvg = 2000L,
    k_min = 3L,
    k_max = NULL,
    batch_frac_threshold_for_k_max = 0.01,
    min_batch_coverage = NULL,
    min_cells_in_domain_per_batch = 5L,
    max_underrepresentation_fold = 10L,
    max_residual_cell_frac = 0.05,
    exclude_cells = NULL,
    seed = 236L,
    output_level = "standard",
    verbose = TRUE
  ) {
    # Resolve serialized S4 classes before dispatch so their load hooks can run.
    if (isS4(object)) methods::getClassDef(class(object))
    standardGeneric("bcmp")
  }
)
methods::setMethod(
  "bcmp",
  signature(object = "ANY"),
  function(
    object,
    batch_key = "batch",
    assay = NULL,
    layer = NULL,
    partition_n_pcs = 30L,
    partition_n_hvg = 2000L,
    k_min = 3L,
    k_max = NULL,
    batch_frac_threshold_for_k_max = 0.01,
    min_batch_coverage = NULL,
    min_cells_in_domain_per_batch = 5L,
    max_underrepresentation_fold = 10L,
    max_residual_cell_frac = 0.05,
    exclude_cells = NULL,
    seed = 236L,
    output_level = "standard",
    verbose = TRUE
  ) {
    stop("`bcmp()` supports only Seurat v5 and SingleCellExperiment objects.", call. = FALSE)
  }
)
