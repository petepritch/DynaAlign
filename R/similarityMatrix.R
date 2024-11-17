# Set number of OpenMP threads
Sys.setenv("OMP_NUM_THREADS" = parallel::detectCores())

#' @importFrom Rcpp sourceCpp
NULL

#' Calculate Sequence Similarity Matrix
#' 
#' @param sequences A character vector containing sequences
#' @return A numeric matrix of pairwise sequence similarities
#' @export
get_similarity_psa <- function(sequences) {
  if(!is.character(sequences)) {
    stop("Input must be a character vector of sequences")
  }
  calculateSimilarityMatrix(sequences)
}

#' @importFrom Rcpp sourceCpp
NULL

#' Calculate MinHash Similarity Matrix
#' 
#' @param sequences A character vector containing sequences
#' @param k Length of k-mers (default: 3)
#' @param num_hash Number of hash functions (default: 100)
#' @return A numeric matrix of pairwise similarities
#' @export
get_similarity_minhash <- function(sequences, k = 2, num_hash = 50) {
  if(!is.character(sequences)) {
    stop("Input must be a character vector of sequences")
  }
  if(k < 1) {
    stop("k must be positive")
  }
  if(num_hash < 1) {
    stop("num_hash must be positive")
  }
  minhash_similarity_matrix(sequences, k, num_hash)
}