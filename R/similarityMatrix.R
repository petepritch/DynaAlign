# Set number of OpenMP threads
Sys.setenv("OMP_NUM_THREADS" = parallel::detectCores())

#' @importFrom Rcpp sourceCpp
NULL

#' Calculate Pairwise Sequence Similarity Matrix
#' 
#' @param sequences A character vector containing peptide sequences.
#' @param matrix_type The substitution matrix to use. Options are "BLOSUM45", 
#'                   "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "BLOSUM100". 
#'                   Default is "BLOSUM62".
#' @param gap_open Gap opening penalty. Default is 10.
#' @param gap_ext Gap extension penalty. Default is 4.
#' 
#' @return A numeric matrix of pairwise sequence similarities, with values 
#'         normalized between 0 and 1.
#' 
#' @details This function performs pairwise sequence alignment using the 
#'          Needleman-Wunsch algorithm with the specified scoring matrix and 
#'          gap penalties. Results are normalized to provide similarity scores 
#'          between 0 and 1.
#' 
#' @examples
#' sequences <- c("HEAGAWGHEE", "PAWHEAE", "HEAGAWGHEE")
#' sim_matrix <- get_similarity_psa(sequences)
#' 
#' @export
get_similarity_psa <- function(sequences, 
                               matrix_type = "BLOSUM62",
                               gap_open = 10, 
                               gap_ext = 4) {
  if(!is.character(sequences)) {
    stop("Input must be a character vector of sequences")
  }
  sim_matrix <- calculateSimilarityMatrix(sequences, 
                                          matrixName = matrix_type,
                                          gapOpen = gap_open, 
                                          gapExt = gap_ext)
  return(sim_matrix)
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