# File: R/plotting.R
#' Plot MinHash Distance Matrix
#'
#' Creates a heatmap visualization of the MinHash distance matrix.
#'
#' @param x A minhash_result object
#' @param cluster logical; whether to cluster the sequences
#' @param ... Additional parameters passed to heatmap
#'
#' @return None (plots to current device)
#' @export
#'
#' @importFrom stats hclust dist
#' @importFrom graphics heatmap
plot.minhash_result <- function(x, cluster = TRUE, ...) {
  heatmap(x$dist_matrix,
          Rowv = if(cluster) as.dendrogram(hclust(dist(x$dist_matrix))) else NA,
          Colv = if(cluster) as.dendrogram(hclust(dist(t(x$dist_matrix)))) else NA,
          main = "Sequence Similarity Heatmap",
          xlab = "Sequence Index",
          ylab = "Sequence Index",
          ...)
}
