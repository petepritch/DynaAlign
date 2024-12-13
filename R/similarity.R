# File: R/similarity.R
#' Calculate Similarity Statistics
#'
#' Computes various similarity statistics for a given similarity matrix.
#'
#' @param X A symmetric similarity matrix
#' @return A list of similarity statistics
#' @export
#' 
#' @importFrom stats median
compute_similarity_stats <- function(X) {
  if (!is.matrix(X)) {
    stop("Input must be a matrix")
  }
  if (!isSymmetric(X)) {
    warning("Input matrix is not symmetric. Results may be unexpected.")
  }
  
  upper_tri_values <- X[upper.tri(X)]
  
  stats <- list(
    mean_similarity = mean(upper_tri_values),
    median_similarity = median(upper_tri_values),
    min_similarity = min(upper_tri_values),
    max_similarity = max(upper_tri_values),
    most_similar_pair = which(X == max(upper_tri_values), 
                              arr.ind = TRUE)[1,],
    least_similar_pair = which(X == min(upper_tri_values), 
                               arr.ind = TRUE)[1,]
  )
  
  class(stats) <- "similarity_stats"
  return(stats)
}