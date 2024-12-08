# File: R/minHash.R
#' Generate k-shingles from a string
#'
#' @param x A character string to be shingled
#' @param k Integer, specifying the length of each shingle
#'
#' @return A character vector containing all k-shingles from the input string
#' @export
#'
#' @examples
#' shingle("ABCDEF", 3) # Returns c("ABC", "BCD", "CDE", "DEF")
shingle <- function(x, k) {
  if (!is.character(x) || length(x) != 1)
    stop("Input 'x' must be a single character string", call. = FALSE)
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k > nchar(x))
    stop(sprintf("'k' must be a positive integer between 1 and %d", nchar(x)), call. = FALSE)
  n <- nchar(x)
  shingles <- vector("character", length = n - k + 1)
  for (i in 1:(n - k + 1)) {
    shingles[i] <- substr(x, i, i + k - 1)
  }
  return(shingles)
}
 
#' Create a vocabulary from sequences
#' 
#' Creates a sorted vocabulary of all unique k-shingles present across all input sequences.
#'
#' @param sequences Character vector of input sequences
#' @param k Integer, specifying the length of each shingle
#'
#' @return Character vector containing all unique k-shingles sorted alphabetically
#' @export
#'
#' @examples 
#' sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
#' create_vocab(sequences, k = 3)
create_vocab <- function(sequences, k) {
  all_shingles <- unique(unlist(lapply(sequences, shingle, k = k)))
  return(sort(all_shingles))
}

#' Create characteristic matrix
#' 
#' Creates a binary matrix where rows represent shingles from the vocabulary
#' and columns represent sequences. A value of 1 indicates the presence of a
#' shingle in a sequence.
#'
#' @param sequences Character vector of input sequences
#' @param vocab Character vector of vocabulary
#' @param k Integer, specifying length of each shingle 
#'
#' @return A binary matrix with dimensions length(vocab) x length(sequences)
#' @export
#'
#' @examples
#' sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
#' vocab <- create_vocab(sequences, k = 3)
#' create_char_matrix(sequences, vocab, k = 3)
create_char_matrix <- function(sequences, vocab, k) {
  seq_shingles <- lapply(sequences, shingle, k = k)
  char_matrix <- sapply(seq_shingles, function(shingles) {
    as.integer(vocab %in% shingles)
  })
  return(char_matrix)
}

#' Generate hash function parameters
#' 
#' Generates parameters for multiple hash functions of the form (ax + b) mod m.
#'
#' @param n_hash Integer number of hash functions to generate
#' @param max_val Integer maximum value for hash function parameters
#'
#' @return List containing two vectors:
#'    \item{a}{Vector of 'a' coefficients for hash functions}
#'    \item{b}{Vector of 'b' coefficients for hash functions}
#' @export
#'
#' @examples create_hash_parameters(n_hash = 10, max_val = 100)
create_hash_parameters <- function(n_hash, max_val) {
  if (n_hash < 1) stop("Number of hash functions must be positive")
  if (max_val < 2) stop("Maximum value must be at least 2")
  
  a_values <- sample(1:max_val, n_hash, replace = TRUE)
  b_values <- sample(0:max_val, n_hash, replace = TRUE)
  return(list(a = a_values, b = b_values))
}

#' Apply a single hash function
#' 
#' Applies a linear hash function of the form (ax + b) mod m to an input value.
#'
#' @param x Integer input value to hash
#' @param a Integer coefficient for linear hash function
#' @param b Integer offset for linear hash function
#' @param m Integer modulus for hash function
#'
#' @return Integer hash value
#' @export
#'
#' @examples
#' apply_hash(x = 5, a = 2, b = 3, m = 100)
apply_hash <- function(x, a, b, m) {  
  return((a * x + b) %% m)
}

#' Compute MinHash signature matrix
#' 
#' Creates signature matrix using multiple hash functions, where each column
#' represents a sequence and each row represents a hash function.
#'
#' @param char_matrix Binary Characteristic matrix
#' @param hash_params List of hash function parameters from create_hash_parameters()
#' @param max_val Integer maximum value for hash functions
#'
#' @return Matrix where each column is a MinHash signature for a sequence
#' @export
#'
#' @examples
#' sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
#' vocab <- create_vocab(sequences, k = 3)
#' char_matrix <- create_char_matrix(sequences, vocab, k = 3)
#' hash_params <- create_hash_parameters(n_hash = 10, max_val = length(vocab))
#' compute_signature_matrix(char_matrix, hash_params, max_val = length(vocab))+
compute_signature_matrix <- function(char_matrix, hash_params, max_val) {
  n_hash <- length(hash_params$a)
  n_docs <- ncol(char_matrix)
  sig_matrix <- matrix(Inf, nrow = n_hash, ncol = n_docs)
  
  for (i in 1:nrow(char_matrix)) {
    hash_values <- mapply(function(a, b) apply_hash(i, a, b, max_val),
                          hash_params$a, hash_params$b)
    
    for (j in 1:n_docs) {
      if (char_matrix[i,j] == 1) {
        sig_matrix[,j] <- pmin(sig_matrix[,j], hash_values)
      }
    }
  }
  
  return(sig_matrix)
}

#' Compute Jaccard distance matrix
#' 
#' Computes pairwise Jaccard distances between sequences using their MinHash signatures.
#'
#' @param sig_matrix Symmetric matrix of Jaccard distances between sequences
#'
#' @return Jaccard distance matrix
#' @export
#'
#' @examples
#' sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
#' k <- 3
#' n_hash <- 10
#' # Create vocabulary and matrices
#' vocab <- create_vocab(sequences, k)
#' char_matrix <- create_char_matrix(sequences, vocab, k)
#' max_val <- length(vocab)
#' hash_params <- create_hash_parameters(n_hash, max_val)
#' sig_matrix <- compute_signature_matrix(char_matrix, hash_params, max_val)
#' # Compute distance matrix
#' compute_distance_matrix(sig_matrix)
compute_distance_matrix <- function(sig_matrix) {
  
  n_docs <- ncol(sig_matrix)
  dist_matrix <- matrix(0, nrow = n_docs, ncol = n_docs)
  
  for (i in 1:n_docs) {
    for (j in i:n_docs) {
      if (i != j) {
        similarity <- mean(sig_matrix[,i] == sig_matrix[,j])
        dist_matrix[i,j] <- 1 - similarity
        dist_matrix[j,i] <- dist_matrix[i,j]  
      }
    }
  }
  
  return(dist_matrix)
}

#' MinHash pipeline for sequence similarity
#' 
#' Complete pipeline for computing Jaccard distances between sequences using MinHash.
#' This function handles all steps from creating shingles to computing the final
#' distance matrix.
#'
#' @param sequences Character vector of input sequences
#' @param k Integer length of shingles
#' @param n_hash Integer number of hash functions to use
#'
#' @return List containing:
#'   \item{vocabulary}{Character vector of all unique k-shingles}
#'   \item{char_matrix}{Binary characteristic matrix}
#'   \item{sig_matrix}{MinHash signature matrix}
#'   \item{dist_matrix}{Jaccard distance matrix}
#' @export
#'
#' @examples
#' sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN", "XXXXXYYYYYYZZ")
#' result <- minhash(sequences, k = 3, n_hash = 100)
#' # View distance matrix
#' print(result$dist_matrix)
minhash <- function(sequences, k, n_hash) {

  vocab <- create_vocab(sequences, k)
  char_matrix <- create_char_matrix(sequences, vocab, k)
  max_val <- length(vocab)
  hash_params <- create_hash_parameters(n_hash, max_val)
  sig_matrix <- compute_signature_matrix(char_matrix, hash_params, max_val)
  dist_matrix <- compute_distance_matrix(sig_matrix)
  
  return(list(
    vocabulary = vocab,
    char_matrix = char_matrix,
    sig_matrix = sig_matrix,
    dist_matrix = dist_matrix
  ))
}
