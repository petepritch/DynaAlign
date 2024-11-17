#' @importFrom Rcpp sourceCpp
NULL

#' Calculate Sequence Similarity Matrix
#' 
#' @param sequences A character vector containing sequences
#' @return A numeric matrix of pairwise sequence similarities
#' @export
get_similarity_matrix <- function(sequences) {
  if(!is.character(sequences)) {
    stop("Input must be a character vector of sequences")
  }
  calculateSimilarityMatrix(sequences)
}