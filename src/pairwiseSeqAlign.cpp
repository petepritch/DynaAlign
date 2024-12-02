#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>

using namespace std;
using namespace Rcpp;

// Amino acid to index mapping (defined once)
const map<char, int> aa_to_index = {
  {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4},
  {'Q', 5}, {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9},
  {'L',10}, {'K',11}, {'M',12}, {'F',13}, {'P',14},
  {'S',15}, {'T',16}, {'W',17}, {'Y',18}, {'V',19}
};

// BLOSUM62 substitution matrix (defined once)
const int blosum62[20][20] = {
  { 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0},  // A
  {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3},  // R
  {-2, 0, 6, 1,-3, 0, 0,-1, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3},  // N
  {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3},  // D
  { 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},  // C
  {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2},  // Q
  {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2},  // E
  { 0,-2,-1,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3},  // G
  {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3},  // H
  {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3},  // I
  {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1},  // L
  {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2},  // K
  {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1},  // M
  {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1},  // F
  {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2},  // P
  { 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2},  // S
  { 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0},  // T
  {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3},  // W
  {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1},  // Y
  { 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4},  // V
};

// [[Rcpp::export]]
double needleman_wunsch_score(const string &sequence1, const string &sequence2,
                              int gapOpen = 10, int gapExt = 4) {
  size_t m = sequence1.size();
  size_t n = sequence2.size();
  
  // Initialize matrices
  vector<vector<int>> M(m + 1, vector<int>(n + 1, 0));
  vector<vector<int>> Ix(m + 1, vector<int>(n + 1, std::numeric_limits<int>::min() / 2));
  vector<vector<int>> Iy(m + 1, vector<int>(n + 1, std::numeric_limits<int>::min() / 2));
  
  // Initialize first row and column
  M[0][0] = 0;
  Ix[0][0] = Iy[0][0] = std::numeric_limits<int>::min() / 2;
  for (size_t i = 1; i <= m; ++i) {
    M[i][0] = std::numeric_limits<int>::min() / 2;
    Ix[i][0] = -gapOpen - (i - 1) * gapExt;
    Iy[i][0] = std::numeric_limits<int>::min() / 2;
  }
  for (size_t j = 1; j <= n; ++j) {
    M[0][j] = std::numeric_limits<int>::min() / 2;
    Ix[0][j] = std::numeric_limits<int>::min() / 2;
    Iy[0][j] = -gapOpen - (j - 1) * gapExt;
  }
  
  // Fill matrices
  for (size_t i = 1; i <= m; ++i) {
    char aa1 = sequence1[i - 1];
    auto it1 = aa_to_index.find(aa1);
    if (it1 == aa_to_index.end()) {
      Rcpp::stop("Invalid amino acid in sequence1: %c", aa1);
    }
    int index1 = it1->second;
    for (size_t j = 1; j <= n; ++j) {
      char aa2 = sequence2[j - 1];
      auto it2 = aa_to_index.find(aa2);
      if (it2 == aa_to_index.end()) {
        Rcpp::stop("Invalid amino acid in sequence2: %c", aa2);
      }
      int index2 = it2->second;
      int score = blosum62[index1][index2];
      
      // Compute Ix[i][j]
      Ix[i][j] = std::max(
        M[i - 1][j] - (gapOpen + gapExt),
        Ix[i - 1][j] - gapExt
      );
      
      // Compute Iy[i][j]
      Iy[i][j] = std::max(
        M[i][j - 1] - (gapOpen + gapExt),
        Iy[i][j - 1] - gapExt
      );
      
      // Compute M[i][j]
      M[i][j] = std::max({
        M[i - 1][j - 1] + score,
        Ix[i - 1][j - 1] + score,
        Iy[i - 1][j - 1] + score
      });
    }
  }
  
  // The optimal alignment score
  int alignment_score = std::max({ M[m][n], Ix[m][n], Iy[m][n] });
  
  return static_cast<double>(alignment_score);
}

// [[Rcpp::export]]
NumericMatrix calculateSimilarityMatrix(CharacterVector sequences,
                                        int gapOpen = 10, int gapExt = 4) {
  size_t n = sequences.length();
  NumericMatrix similarityMatrix(n, n);
  NumericVector selfScores(n);
  
  // Compute self-alignment scores
  for (size_t i = 0; i < n; ++i) {
    string seq = as<string>(sequences[i]);
    double score = needleman_wunsch_score(seq, seq, gapOpen, gapExt);
    selfScores[i] = score;
  }
  
  // Adjust scores if negative
  double minSelfScore = *std::min_element(selfScores.begin(), selfScores.end());
  double adjustValue = (minSelfScore < 0) ? (-minSelfScore + 1) : 0;
  
  for (size_t i = 0; i < n; ++i) {
    selfScores[i] += adjustValue;
  }
  
  // Calculate pairwise similarities
  for (size_t i = 0; i < n; ++i) {
    string seq1 = as<string>(sequences[i]);
    for (size_t j = i; j < n; ++j) {
      if (i == j) {
        similarityMatrix(i, j) = 1.0;  // Diagonal elements
        continue;
      }
      
      string seq2 = as<string>(sequences[j]);
      
      // Get alignment score
      double score = needleman_wunsch_score(seq1, seq2, gapOpen, gapExt);
      score += adjustValue;
      
      // Normalize the score
      double normScore = score / sqrt(selfScores[i] * selfScores[j]);
      
      // Ensure the similarity is between 0 and 1
      normScore = std::max(0.0, std::min(1.0, normScore));
      
      // Fill both upper and lower triangles
      similarityMatrix(i, j) = normScore;
      similarityMatrix(j, i) = normScore;
    }
  }
  
  // Add dimension names
  // Create numeric labels 1,2,3...
  CharacterVector labels(n);
  for(size_t i = 0; i < n; ++i) {
    labels[i] = std::to_string(i + 1);
  }
  
  // Add dimension names using numeric labels
  similarityMatrix.attr("dimnames") = List::create(labels, labels);
  
  return similarityMatrix;
}