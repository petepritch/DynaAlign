# Test shingle function
test_that("shingle function works correctly", {
  # Test basic functionality
  expect_equal(shingle("ABCDEF", 3), c("ABC", "BCD", "CDE", "DEF"))
  
  # Test error handling
  expect_error(shingle(123, 3), "Input 'x' must be a single character string")
  expect_error(shingle("ABCDEF", 0), "'k' must be a positive integer between 1 and 6")
  expect_error(shingle("ABCDEF", 7), "'k' must be a positive integer between 1 and 6")
  
  # Test edge cases
  expect_equal(shingle("AB", 2), "AB")
  expect_length(shingle("ABCDEF", 1), 6)
})

# Test create_vocab function
test_that("create_vocab function works correctly", {
  # Test basic functionality
  sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
  vocab <- create_vocab(sequences, k = 3)
  
  # Check that vocab is sorted
  expect_true(all(vocab == sort(vocab)))
  
  # Check unique shingles
  expect_true(length(vocab) == length(unique(vocab)))
  
  # Verify shingles are correct length
  expect_true(all(nchar(vocab) == 3))
})

# Test create_char_matrix function
test_that("create_char_matrix function works correctly", {
  sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
  vocab <- create_vocab(sequences, k = 3)
  char_matrix <- create_char_matrix(sequences, vocab, k = 3)
  
  # Check matrix dimensions
  expect_equal(nrow(char_matrix), length(vocab))
  expect_equal(ncol(char_matrix), length(sequences))
  
  # Check binary matrix properties
  expect_true(all(char_matrix %in% c(0, 1)))
})

# Test create_hash_parameters function
test_that("create_hash_parameters function works correctly", {
  # Test generation of hash parameters
  n_hash <- 10
  max_val <- 100
  hash_params <- create_hash_parameters(n_hash, max_val)
  
  # Check correct number of parameters
  expect_length(hash_params$a, n_hash)
  expect_length(hash_params$b, n_hash)
  
  # Check parameter ranges
  expect_true(all(hash_params$a > 0 & hash_params$a <= max_val))
  expect_true(all(hash_params$b >= 0 & hash_params$b <= max_val))
})

# Test apply_hash function
test_that("apply_hash function works correctly", {
  # Test basic hash application
  result <- apply_hash(x = 5, a = 2, b = 3, m = 100)
  
  # Check result is within modulus
  expect_true(result >= 0 && result < 100)
  
  # Test deterministic behavior
  expect_equal(apply_hash(5, 2, 3, 100), apply_hash(5, 2, 3, 100))
})

# Test compute_signature_matrix function
test_that("compute_signature_matrix function works correctly", {
  sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN")
  vocab <- create_vocab(sequences, k = 3)
  char_matrix <- create_char_matrix(sequences, vocab, k = 3)
  hash_params <- create_hash_parameters(n_hash = 10, max_val = length(vocab))
  
  sig_matrix <- compute_signature_matrix(char_matrix, hash_params, max_val = length(vocab))
  
  # Check matrix dimensions
  expect_equal(nrow(sig_matrix), length(hash_params$a))
  expect_equal(ncol(sig_matrix), length(sequences))
  
  # Check signature values are numeric
  expect_true(is.numeric(sig_matrix))
})

# Test compute_distance_matrix function
test_that("compute_distance_matrix function works correctly", {
  # Create a mock signature matrix
  sig_matrix <- matrix(c(
    1, 2, 3,
    1, 2, 4,
    2, 3, 5
  ), nrow = 3, ncol = 3)
  
  dist_matrix <- compute_distance_matrix(sig_matrix)
  
  # Check matrix properties
  expect_true(isSymmetric(dist_matrix))
  expect_equal(diag(dist_matrix), rep(0, ncol(dist_matrix)))
  expect_true(all(dist_matrix >= 0 & dist_matrix <= 1))
})

# Test minhash pipeline
test_that("minhash pipeline works end-to-end", {
  sequences <- c("ACDEGHHIKLLL", "ACDEGHHIKLMN", "XXXXXYYYYYYZZ")
  result <- minhash(sequences, k = 3, n_hash = 100)
  
  # Check returned list components
  expect_true(all(c("vocabulary", "char_matrix", "sig_matrix", "dist_matrix") %in% names(result)))
  
  # Verify matrix dimensions
  expect_equal(nrow(result$char_matrix), length(result$vocabulary))
  expect_equal(ncol(result$char_matrix), length(sequences))
  expect_equal(nrow(result$sig_matrix), 100)
  expect_equal(ncol(result$sig_matrix), length(sequences))
  expect_equal(dim(result$dist_matrix), c(length(sequences), length(sequences)))
})