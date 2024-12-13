# File: R/plotting.R
#' Plot Similarity Matrix
#'
#' Creates a heatmap visualization of a similarity matrix.
#'
#' @param X A similarity matrix (symmetric matrix of pairwise similarities)
#' @param cluster logical; whether to cluster the sequences
#' @param ... Additional parameters passed to heatmap
#'
#' @return None (plots to current device)
#' @export
#'
#' @importFrom stats hclust dist as.dendrogram heatmap
plot_similarity_matrix <- function(X, cluster = TRUE, ...) {
  if (!is.matrix(X)) {
    stop("Input must be a matrix")
  }
  if (!isSymmetric(X)) {
    warning("Input matrix is not symmetric. Results may be unexpected.")
  }
  
  heatmap(X,
          Rowv = if(cluster) as.dendrogram(hclust(dist(X))) else NA,
          Colv = if(cluster) as.dendrogram(hclust(dist(t(X)))) else NA,
          main = "Similarity Matrix Heatmap",
          xlab = "Sequence/Item Index",
          ylab = "Sequence/Item Index",
          ...)
}