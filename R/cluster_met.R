#' Cluster cells in a methylTracer object
#'
#' @description
#' Perform clustering of cells based on HDF5-backed region-level methylation
#' values stored in a \code{methylTracer} object. The function implements a
#' standard workflow:
#' \enumerate{
#'   \item Select highly variable features.
#'   \item Perform PCA on the selected features.
#'   \item Run k-means clustering in PCA space and write cluster labels to
#'         \code{/obs/<cluster_colname>}.
#'   \item Optionally compute a 2D UMAP embedding and save it under
#'         \code{/uns/umap_coords}.
#'   \item Optionally build a k-nearest neighbor graph from the UMAP embedding
#'         and perform Leiden community detection, writing labels to
#'         \code{/obs/<leiden_colname>}.
#' }
#'
#' The function can either require the input matrix to be free of missing
#' values, or it can perform a simple feature-wise mean imputation internally.
#' For single-cell methylation data, it is generally recommended to run
#' \code{compute_qc_value()}, \code{filter_obs_var()} and \code{impute_met_obj()}
#' before clustering, and use \code{na_method = "error"} here as a safeguard.
#'
#' @param met A \code{methylTracer} object with an on-disk methylation matrix
#'   stored under dataset \code{"X"}.
#' @param n_clusters Integer scalar; number of k-means clusters in PCA space.
#' @param top_features Integer scalar; number of features (rows of \code{X})
#'   with the highest variance (across cells) to use for PCA and clustering.
#' @param n_pcs Integer scalar; number of principal components to retain from
#'   \code{\link[stats]{prcomp}} and use for k-means, UMAP and Leiden.
#' @param cluster_colname Character scalar; name of the column written to
#'   \code{/obs/<cluster_colname>} with k-means cluster labels. Defaults to
#'   \code{"cluster"}.
#' @param seed Integer scalar; random seed used for PCA initialization,
#'   k-means, UMAP and Leiden (via underlying stochastic routines).
#' @param save_pca Logical; if \code{TRUE}, PCA coordinates are written to
#'   \code{/uns/pca_coords} (or the dataset specified by \code{pca_dataset}).
#' @param pca_dataset Character scalar; dataset path under \code{/uns} where
#'   PCA coordinates are stored. Defaults to \code{"uns/pca_coords"}.
#' @param run_umap Logical; if \code{TRUE}, compute a 2D UMAP embedding from
#'   the retained PCs.
#' @param umap_n_neighbors Integer scalar; number of neighbors for UMAP.
#'   Larger values emphasize global structure at the cost of local detail.
#' @param umap_min_dist Numeric scalar; UMAP \code{min_dist} parameter controlling
#'   how tightly points are allowed to pack in the embedding space.
#' @param umap_dataset Character scalar; dataset path under \code{/uns} where
#'   UMAP coordinates are stored. Defaults to \code{"uns/umap_coords"}.
#' @param run_leiden Logical; if \code{TRUE}, construct a k-nearest neighbor
#'   graph from the UMAP nearest neighbor output and perform Leiden clustering.
#' @param leiden_resolution Numeric scalar; resolution parameter passed to
#'   \code{\link[igraph]{cluster_leiden}}. Higher values typically yield more
#'   fine-grained clusters.
#' @param leiden_colname Character scalar; name of the column written to
#'   \code{/obs/<leiden_colname>} with Leiden cluster labels. Defaults to
#'   \code{"leiden"}.
#' @param na_method Character scalar; how to handle missing values in the
#'   selected feature matrix. One of:
#'   \itemize{
#'     \item \code{"error"}: stop with an error if any NA is present.
#'     \item \code{"feature_mean"}: drop features with NA fraction greater
#'           than \code{max_na_frac} and mean-impute remaining NA values
#'           per feature (column).
#'   }
#' @param max_na_frac Numeric scalar in \eqn{[0,1]}; when
#'   \code{na_method = "feature_mean"}, features with a fraction of NA values
#'   greater than \code{max_na_frac} are excluded before PCA. Defaults to
#'   \code{0.5}.
#'
#' @details
#' The function operates directly on the HDF5-backed matrix \code{X} via
#' \pkg{HDF5Array} and \pkg{DelayedMatrixStats}. Variance-based feature
#' selection is used to restrict computation to \code{top_features} rows.
#' PCA is then performed on the (possibly imputed) in-memory matrix of
#' shape \code{cells x features}.
#'
#' UMAP is computed using \pkg{uwot} on the retained PCs. When
#' \code{run_leiden = TRUE}, the k-nearest neighbor information returned by
#' UMAP (\code{ret_nn = TRUE}) is used to build an undirected graph (via
#' \pkg{igraph}), on which Leiden community detection is performed.
#'
#' K-means cluster labels (\code{"C1", "C2", ...}) are written to
#' \code{/obs/<cluster_colname>}, and Leiden cluster labels (\code{"L1", "L2", ...})
#' are written to \code{/obs/<leiden_colname>} in the underlying HDF5 file.
#' PCA and UMAP coordinates are stored under \code{/uns} if requested.
#'
#' @return
#' A named list with components:
#' \itemize{
#'   \item \code{cluster}: character vector of k-means cluster labels
#'         (length = number of cells).
#'   \item \code{pca}: data.frame of PCA coordinates (cells x PCs).
#'   \item \code{umap}: data.frame of UMAP coordinates (cells x 2), or
#'         \code{NULL} if \code{run_umap = FALSE}.
#'   \item \code{leiden}: character vector of Leiden cluster labels, or
#'         \code{NULL} if \code{run_leiden = FALSE}.
#' }
#'
#' Side effects include writing cluster labels and embeddings into the HDF5
#' file associated with \code{met}.
#'
#' @seealso
#' \code{\link{compute_qc_value}}, \code{\link{filter_obs_var}},
#' \code{\link{impute_met_obj}}, \code{\link[uwot]{umap}},
#' \code{\link[igraph]{cluster_leiden}}
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedMatrixStats rowVars
#' @importFrom stats prcomp kmeans
#' @importFrom rhdf5 h5write h5delete
#' @importFrom uwot umap
#' @importFrom igraph graph_from_data_frame cluster_leiden
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'met' is a methylTracer object that has been QC-filtered
#' # and imputed (e.g. via compute_qc_value(), filter_obs_var(), impute_met_obj()).
#'
#' clust_res <- cluster_met_cells(
#'   met,
#'   n_clusters       = 10,
#'   top_features     = 2000,
#'   n_pcs            = 30,
#'   run_umap         = TRUE,
#'   umap_n_neighbors = 30,
#'   umap_min_dist    = 0.3,
#'   run_leiden       = TRUE,
#'   leiden_resolution = 1.0,
#'   na_method        = "error"
#' )
#'
#' table(clust_res$cluster)
#' table(clust_res$leiden)
#' }
cluster_met_cells <- function(
    met,
    n_clusters       = 2,
    top_features     = 2000,
    n_pcs            = 20,
    cluster_colname  = "cluster",
    seed             = 1,
    save_pca         = TRUE,
    pca_dataset      = "uns/pca_coords",
    run_umap         = TRUE,
    umap_n_neighbors = 3,
    umap_min_dist    = 0.1,
    umap_dataset     = "uns/umap_coords",
    run_leiden       = TRUE,
    leiden_resolution = 0.5,
    leiden_colname    = "leiden",
    ## how to handle NA values
    na_method        = c("error", "feature_mean"),
    max_na_frac      = 0.5  # drop features with > 50% NA
) {
  na_method <- match.arg(na_method)
  
  h5file     <- met@seed@filepath
  cell_names <- as.vector(met@sample_name)
  
  ## 1) HDF5-backed matrix
  X <- HDF5Array::HDF5Array(h5file, name = "X")
  
  ## 2) Select highly variable features (allowing NA via na.rm = TRUE)
  vars <- DelayedMatrixStats::rowVars(X, na.rm = TRUE)
  o    <- order(vars, decreasing = TRUE)
  k    <- min(top_features, length(o))
  sel_idx <- o[seq_len(k)]
  
  X_sel <- X[sel_idx, ]
  mat   <- t(as.matrix(X_sel))  # cells x features
  
  ## 2.5) Handle NA values if present
  if (anyNA(mat)) {
    if (na_method == "error") {
      stop(
        "The selected feature matrix for clustering contains NA values. ",
        "Please run 'impute_met_obj()' first to remove NA, or call ",
        "cluster_met_cells() with na_method = 'feature_mean'."
      )
    }
    
    # na_method == "feature_mean"
    # (1) Drop features with too high NA fraction
    na_frac <- colMeans(is.na(mat))
    keep_feature <- na_frac <= max_na_frac
    
    if (!any(keep_feature)) {
      stop(
        "All selected features have NA fraction > max_na_frac (",
        max_na_frac, "). Try a larger max_na_frac or impute first."
      )
    }
    
    if (any(!keep_feature)) {
      mat <- mat[, keep_feature, drop = FALSE]
      if (ncol(mat) < 2L) {
        stop("Too few features remain after NA filtering to run PCA.")
      }
    }
    
    # (2) Mean-impute remaining NA values per feature (column)
    col_means <- colMeans(mat, na.rm = TRUE)
    for (j in seq_len(ncol(mat))) {
      idx_na <- is.na(mat[, j])
      if (any(idx_na)) {
        mat[idx_na, j] <- col_means[j]
      }
    }
  }
  
  ## 3) PCA on selected (and possibly imputed) features
  set.seed(seed)
  pca <- stats::prcomp(mat, center = TRUE, scale. = FALSE)
  n_use_pcs <- min(n_pcs, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_use_pcs), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(n_use_pcs))
  rownames(pcs) <- cell_names
  
  ## 4) k-means clustering on PCA space
  km <- stats::kmeans(pcs, centers = n_clusters, nstart = 10)
  clusters <- paste0("C", km$cluster)
  
  ## 5) Write k-means clusters to /obs/<cluster_colname>
  cluster_path <- paste0("obs/", cluster_colname)
  try(rhdf5::h5delete(h5file, cluster_path), silent = TRUE)
  rhdf5::h5write(
    obj  = clusters,
    file = h5file,
    name = cluster_path
  )
  
  ## 6) Optionally save PCA coordinates
  if (save_pca) {
    try(rhdf5::h5delete(h5file, pca_dataset), silent = TRUE)
    rhdf5::h5write(
      obj  = pcs,
      file = h5file,
      name = pca_dataset
    )
  }
  pca_df <- as.data.frame(pcs, check.names = FALSE)
  
  ## 7) Optionally run UMAP
  umap_df <- NULL
  umap_nn <- NULL
  if (run_umap || run_leiden) {
    set.seed(seed)
    umap_res <- uwot::umap(
      pcs,
      n_neighbors = umap_n_neighbors,
      min_dist    = umap_min_dist,
      n_components = 2,
      ret_nn      = run_leiden,
      verbose     = FALSE
    )
    
    if (run_leiden) {
      emb <- umap_res$embedding
      nn  <- umap_res$nn
    } else {
      emb <- umap_res
      nn  <- NULL
    }
    
    colnames(emb) <- c("UMAP1", "UMAP2")
    rownames(emb) <- cell_names
    
    if (run_umap) {
      try(rhdf5::h5delete(h5file, umap_dataset), silent = TRUE)
      rhdf5::h5write(
        obj  = emb,
        file = h5file,
        name = umap_dataset
      )
      umap_df <- as.data.frame(emb, check.names = FALSE)
    }
    
    umap_nn <- nn
  }
  
  ## 8) Optionally run Leiden clustering on the UMAP nearest neighbor graph
  leiden_clusters <- NULL
  if (run_leiden) {
    if (is.null(umap_nn)) {
      stop("Leiden clustering requires UMAP nearest neighbor information. ",
           "Please set run_umap = TRUE or ensure ret_nn is available.")
    }
    
    # uwot::umap with ret_nn = TRUE typically returns:
    # nn$idx: integer matrix (cells x neighbors)
    # nn$dist: numeric matrix (cells x neighbors)
    if (!is.null(umap_nn$idx) && !is.null(umap_nn$dist)) {
      nn_idx  <- umap_nn$idx
      nn_dist <- umap_nn$dist
    } else if (!is.null(umap_nn$euclidean$idx) &&
               !is.null(umap_nn$euclidean$dist)) {
      # compatibility for older uwot versions that nest under $euclidean
      nn_idx  <- umap_nn$euclidean$idx
      nn_dist <- umap_nn$euclidean$dist
    } else {
      stop("Unsupported structure of UMAP nearest neighbor output.")
    }
    
    n_cells     <- nrow(nn_idx)
    n_neighbors <- ncol(nn_idx)
    
    from_idx <- rep(seq_len(n_cells), times = n_neighbors)
    to_idx   <- as.vector(nn_idx)
    weight   <- as.vector(nn_dist)
    
    # remove self-loops
    mask <- from_idx != to_idx
    from_idx <- from_idx[mask]
    to_idx   <- to_idx[mask]
    weight   <- weight[mask]
    
    from_cells <- cell_names[from_idx]
    to_cells   <- cell_names[to_idx]
    
    edge_df <- data.frame(
      from   = from_cells,
      to     = to_cells,
      weight = weight,
      stringsAsFactors = FALSE
    )
    
    g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
    cl_obj <- igraph::cluster_leiden(
      g,
      resolution_parameter = leiden_resolution
    )
    
    membership <- cl_obj$membership
    names(membership) <- cl_obj$names
    
    # order Leiden labels in the same order as cell_names
    leiden_raw <- membership[cell_names]
    leiden_clusters <- paste0("L", leiden_raw)
    
    # write Leiden clusters to /obs/<leiden_colname>
    leiden_path <- paste0("obs/", leiden_colname)
    try(rhdf5::h5delete(h5file, leiden_path), silent = TRUE)
    rhdf5::h5write(
      obj  = leiden_clusters,
      file = h5file,
      name = leiden_path
    )
  }
  
  list(
    cluster = clusters,
    pca     = pca_df,
    umap    = umap_df,
    leiden  = leiden_clusters
  )
}
