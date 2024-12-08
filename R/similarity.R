# File: R/similarity.R
#' Calculate Similarity Statistics
#'
#' Computes various similarity statistics for the sequence set.
#'
#' @param result A minhash_result object
#' @return A list of similarity statistics
#' @export
#' 
#' @importFrom stats median
compute_similarity_stats <- function(result) {
  dist_matrix <- result$dist_matrix
  
  stats <- list(
    mean_similarity = mean(1 - dist_matrix[upper.tri(dist_matrix)]),
    median_similarity = median(1 - dist_matrix[upper.tri(dist_matrix)]),
    min_similarity = min(1 - dist_matrix[upper.tri(dist_matrix)]),
    max_similarity = max(1 - dist_matrix[upper.tri(dist_matrix)]),
    most_similar_pair = which(dist_matrix == min(dist_matrix[dist_matrix > 0]), 
                              arr.ind = TRUE)[1,],
    least_similar_pair = which(dist_matrix == max(dist_matrix), 
                               arr.ind = TRUE)[1,]
  )
  
  class(stats) <- "minhash_stats"
  return(stats)
}