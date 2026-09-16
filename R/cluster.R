BCMP_ALL_SINGLETON_CLUSTER_ERROR <-
  "All clusters are singletons after initial clustering; BCMP cannot proceed."

with_bcmp_local_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit(
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    },
    add = TRUE
  )
  set.seed(as.integer(seed))
  force(code)
}

edge_adjacency <- function(n_cells, edges) {
  adjacency <- vector("list", n_cells)
  if (!nrow(edges)) {
    return(adjacency)
  }
  for (edge_index in seq_len(nrow(edges))) {
    left <- edges[edge_index, 1L] + 1L
    right <- edges[edge_index, 2L] + 1L
    adjacency[[left]] <- c(adjacency[[left]], right)
    adjacency[[right]] <- c(adjacency[[right]], left)
  }
  adjacency
}

group_singletons_from_edges <- function(labels, n_cells, edges, seed = 236L) {
  labels <- as.integer(labels)
  label_counts <- table(labels)
  singleton_labels <- as.integer(names(label_counts)[label_counts == 1L])
  if (!length(singleton_labels)) {
    return(as.integer(match(labels, sort(unique(labels))) - 1L))
  }

  cluster_labels <- as.integer(names(label_counts)[label_counts != 1L])
  if (!length(cluster_labels)) {
    stop(BCMP_ALL_SINGLETON_CLUSTER_ERROR, call. = FALSE)
  }

  adjacency <- edge_adjacency(n_cells, edges)
  cluster_cells <- lapply(
    cluster_labels,
    function(label) which(labels == label)
  )
  names(cluster_cells) <- as.character(cluster_labels)

  grouped <- with_bcmp_local_seed(seed, {
    for (singleton_label in singleton_labels) {
      singleton_cell <- which(labels == singleton_label)
      scores <- vapply(
        cluster_labels,
        function(cluster_label) {
          cells <- cluster_cells[[as.character(cluster_label)]]
          sum(adjacency[[singleton_cell]] %in% cells) / length(cells)
        },
        numeric(1)
      )
      candidates <- cluster_labels[scores == max(scores)]
      chosen <- candidates[sample.int(length(candidates), size = 1L)]
      labels[singleton_cell] <- chosen
      cluster_cells[[as.character(chosen)]] <-
        c(cluster_cells[[as.character(chosen)]], singleton_cell)
    }
    labels
  })
  as.integer(match(grouped, sort(unique(grouped))) - 1L)
}

#' Cluster SNN upper-triangle edges with Seurat-compatible Leiden settings.
#'
#' @param n_cells Number of cells.
#' @param edges Integer, zero-based upper-triangle edge matrix.
#' @param seed Local clustering and singleton-grouping seed.
#'
#' @return Zero-based integer domain labels.
#' @noRd
cluster_labels_from_edge_matrix <- function(n_cells, edges, seed = 236L) {
  if (!is.numeric(n_cells) || length(n_cells) != 1L || is.na(n_cells) ||
      !is.finite(n_cells) || n_cells < 1L || n_cells != as.integer(n_cells)) {
    stop("'n_cells' must be one finite positive integer.", call. = FALSE)
  }
  n_cells <- as.integer(n_cells)
  if (!is.matrix(edges) || ncol(edges) != 2L || !is.numeric(edges) ||
      any(!is.finite(edges)) || any(edges != floor(edges))) {
    stop("'edges' must be a finite integer matrix with two columns.", call. = FALSE)
  }
  edges <- matrix(as.integer(edges), ncol = 2L)
  if (nrow(edges) && (any(edges < 0L) || any(edges >= n_cells) ||
      any(edges[, 1L] >= edges[, 2L]))) {
    stop("'edges' must contain valid zero-based upper-triangle pairs.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != as.integer(seed)) {
    stop("'seed' must be one finite integer.", call. = FALSE)
  }

  graph <- igraph::make_empty_graph(n = n_cells, directed = FALSE)
  if (nrow(edges)) {
    graph <- igraph::add_edges(graph, as.vector(t(edges + 1L)))
  }
  membership <- leidenbase::leiden_find_partition(
    graph,
    partition_type = "RBConfigurationVertexPartition",
    edge_weights = NULL,
    resolution_parameter = 0.1,
    seed = max(1L, as.integer(seed)),
    num_iter = 10L
  )$membership
  group_singletons_from_edges(
    labels = membership,
    n_cells = n_cells,
    edges = edges,
    seed = seed
  )
}
