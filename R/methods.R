#' Calculate Similarity Using MinHash
#'
#' This function calculates the similarity of sequences using the MinHash algorithm.
#' It provides a method for quickly estimating similarity between large sets of sequences.
#'
#' @param sequences A character vector of sequences to be compared.
#' @param k An integer specifying the length of each shingle. Default is 4.
#' @param n_hash An integer specifying the number of hash functions to use. Default is 50.
#' 
#' @return A numeric value representing the similarity score between the input sequences.
#' @export
#' 
#' @examples
#' sequences <- c("ATGC", "GTCA")
#' similarityMH(sequences)
similarity.MH <- function(sequences, k = 4L, n_hash = 50L) {
  if (!is.character(sequences) || length(sequences) < 2) {
    stop("Input 'sequences' must be a character vector of at least two sequences.")
  }
  
  if (!is.integer(k) || length(k) != 1 || k <= 0) {
    stop("'k' must be a positive integer.")
  }
  
  if (!is.integer(n_hash) || length(n_hash) != 1 || n_hash <= 0) {
    stop("'n_hash' must be a positive integer.")
  }
  
  .Call('_DynaAlign_similarityMH', PACKAGE = 'DynaAlign', sequences, k, n_hash)
}

#' Calculate Similarity Using Needleman-Wunsch
#'
#' This function calculates the similarity of sequences using the Needleman-Wunsch algorithm.
#' It allows for flexible scoring using a substitution matrix and gap penalties.
#'
#' @param sequences A character vector of sequences to be aligned.
#' @param matrixName A string specifying the substitution matrix to be used. Default is "BLOSUM62".
#' @param gapOpen An integer specifying the gap opening penalty. Default is 10.
#' @param gapExt An integer specifying the gap extension penalty. Default is 4.
#' 
#' @return A numeric value representing the similarity score between the aligned sequences.
#' @export
#' 
#' @examples
#' sequences <- c("ATGC", "GTCA")
#' similarityNW(sequences)
similarity.NW <- function(sequences, matrixName = "BLOSUM62", gapOpen = 10L, gapExt = 4L) {
  if (!is.character(sequences) || length(sequences) != 2) {
    stop("Input 'sequences' must be a character vector of exactly two sequences.")
  }
  
  if (!is.character(matrixName) || length(matrixName) != 1) {
    stop("'matrixName' must be a single character string specifying a substitution matrix.")
  }
  
  if (!is.integer(gapOpen) || length(gapOpen) != 1 || gapOpen < 0) {
    stop("'gapOpen' must be a non-negative integer.")
  }
  
  if (!is.integer(gapExt) || length(gapExt) != 1 || gapExt < 0) {
    stop("'gapExt' must be a non-negative integer.")
  }
  
  .Call('_DynaAlign_similarityNW', PACKAGE = 'DynaAlign', sequences, matrixName, gapOpen, gapExt)
}
