% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minHash.R
\name{compute_signature_matrix}
\alias{compute_signature_matrix}
\title{Compute MinHash signature matrix}
\usage{
compute_signature_matrix(char_matrix, hash_params, max_val)
}
\arguments{
\item{char_matrix}{Binary Characteristic matrix}

\item{hash_params}{List of hash function parameters from create_hash_parameters()}

\item{max_val}{Integer maximum value for hash functions}
}
\value{
Matrix where each column is a MinHash signature for a sequence
}
\description{
Creates signature matrix using multiple hash functions, where each column
represents a sequence and each row represents a hash function.
}
\examples{
sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
vocab <- create_vocab(sequences, k = 3)
char_matrix <- create_char_matrix(sequences, vocab, k = 3)
hash_params <- create_hash_parameters(n_hash = 10, max_val = length(vocab))
compute_signature_matrix(char_matrix, hash_params, max_val = length(vocab))+
}
