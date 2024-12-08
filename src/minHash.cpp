#include <Rcpp.h>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>

// Add OpenMP if available
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

// MurmurHash3 implementation
uint32_t murmur3_32(const char* key, size_t len, uint32_t seed) {
  static const uint32_t c1 = 0xcc9e2d51;
  static const uint32_t c2 = 0x1b873593;
  static const uint32_t r1 = 15;
  static const uint32_t r2 = 13;
  static const uint32_t m = 5;
  static const uint32_t n = 0xe6546b64;
  
  uint32_t hash = seed;
  
  const int nblocks = len / 4;
  const uint32_t* blocks = (const uint32_t*)(key);
  
  for(int i = 0; i < nblocks; i++) {
    uint32_t k = blocks[i];
    k *= c1;
    k = (k << r1) | (k >> (32 - r1));
    k *= c2;
    hash ^= k;
    hash = ((hash << r2) | (hash >> (32 - r2))) * m + n;
  }
  
  const uint8_t* tail = (const uint8_t*)(key + nblocks * 4);
  uint32_t k1 = 0;
  
  switch(len & 3) {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
    k1 *= c1;
    k1 = (k1 << r1) | (k1 >> (32 - r1));
    k1 *= c2;
    hash ^= k1;
  }
  
  hash ^= len;
  hash ^= (hash >> 16);
  hash *= 0x85ebca6b;
  hash ^= (hash >> 13);
  hash *= 0xc2b2ae35;
  hash ^= (hash >> 16);
  
  return hash;
}

// Hash family class for multiple hash functions
class HashFamily {
private:
  vector<uint32_t> seeds;
  
public:
  HashFamily(int num_hash) {
    seeds.resize(num_hash);
    for(int i = 0; i < num_hash; ++i) {
      seeds[i] = i + 1;  // Simple seed generation
    }
  }
  
  uint32_t hash(const string& s, int index) {
    return murmur3_32(s.c_str(), s.length(), seeds[index]);
  }
};

// Generate k-mers from a sequence
vector<string> generate_kmers(const string& seq, int k) {
  vector<string> kmers;
  if (seq.length() >= static_cast<std::size_t>(k)) {
    for(size_t i = 0; i <= seq.length() - k; ++i) {
      kmers.push_back(seq.substr(i, k));
    }
  }
  return kmers;
}

// MinHash Similarity Matrix Computation
//' Compute MinHash Similarity Matrix
 //' 
 //' This function computes a similarity matrix using the MinHash technique
 //' 
 //' @param sequences A character vector of input sequences
 //' @param k The length of k-mers to use (default: 4)
 //' @param n_hash Number of hash functions to use (default: 50)
 //' @return A numeric matrix of pairwise similarities
 //' @export
 // [[Rcpp::export]]
 NumericMatrix similarityMH(CharacterVector sequences, int k = 4, int n_hash = 50) {
   // Input validation
   if (k <= 0) {
     Rcpp::stop("'k' must be a positive integer");
   }
   if (n_hash <= 0) {
     Rcpp::stop("Number of hash functions must be positive");
   }
   
   size_t n = sequences.length();
   NumericMatrix similarityMatrix(n, n);
   
   // Initialize hash family
   HashFamily hash_family(n_hash);
   
   // Store signatures for each sequence
   vector<vector<uint32_t>> signatures(n, vector<uint32_t>(n_hash, UINT32_MAX));
   
   // Parallel processing of signature generation
#ifdef _OPENMP
#pragma omp parallel for
#endif
   for(size_t i = 0; i < n; ++i) {
     string seq = as<string>(sequences[i]);
     vector<string> kmers = generate_kmers(seq, k);
     
     // For each k-mer, update signature
     for(const string& kmer : kmers) {
       for(int h = 0; h < n_hash; ++h) {
         uint32_t hash_value = hash_family.hash(kmer, h);
         signatures[i][h] = min(signatures[i][h], hash_value);
       }
     }
   }
   
   // Calculate similarities
   for(size_t i = 0; i < n; ++i) {
     similarityMatrix(i,i) = 1.0;  // Diagonal elements
     
     // Parallel processing of similarity computation
#ifdef _OPENMP
#pragma omp parallel for
#endif
     for(size_t j = i+1; j < n; ++j) {
       int matches = 0;
       for(int h = 0; h < n_hash; ++h) {
         if(signatures[i][h] == signatures[j][h]) {
           ++matches;
         }
       }
       double similarity = static_cast<double>(matches) / n_hash;
       similarityMatrix(i,j) = similarity;
       similarityMatrix(j,i) = similarity;
     }
   }
   
   // Add dimension names (1,2,3...)
   CharacterVector labels(n);
   for(size_t i = 0; i < n; ++i) {
     labels[i] = to_string(i + 1);
   }
   similarityMatrix.attr("dimnames") = List::create(labels, labels);
   
   return similarityMatrix;
 }

// Additional helper functions can be added here if needed

// Rcpp module registration
#include <Rcpp.h>

RCPP_MODULE(minhash_module) {
  Rcpp:function("similarityMH", &similarityMH, 
           "Compute MinHash similarity matrix for input sequences\n"
           "@param sequences Character vector of input sequences\n"
           "@param k Length of k-mers to use\n"
           "@param n_hash Number of hash functions\n"
           "@return Numeric matrix of pairwise sequence similarities");
}