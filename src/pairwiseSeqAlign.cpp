#include <Rcpp.h>
#include <string>
#include <vector>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
double needleman_wunsch_score(const string &sequence1, const string &sequence2,
                              int matchScore=1, int mismatchPenalty=-2, int gapPenalty=-2) {  
  size_t m = sequence1.size();
  size_t n = sequence2.size();  
  
  // Initialize the dynamic programming matrix  
  vector<vector<int>> dpMatrix(m + 1, vector<int>(n + 1, 0));  
  
  // Initialize the first row and column with gap penalties  
  for (size_t i = 1; i <= m; ++i) {  
    dpMatrix[i][0] = dpMatrix[i - 1][0] + gapPenalty;  
  }  
  
  for (size_t j = 1; j <= n; ++j) {  
    dpMatrix[0][j] = dpMatrix[0][j - 1] + gapPenalty;  
  }  
  
  // Fill the dynamic programming matrix  
  for (size_t i = 1; i <= m; ++i) {  
    for (size_t j = 1; j <= n; ++j) {  
      int match = dpMatrix[i - 1][j - 1] + 
        (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchPenalty);  
      int gap1 = dpMatrix[i - 1][j] + gapPenalty;  
      int gap2 = dpMatrix[i][j - 1] + gapPenalty;  
      
      dpMatrix[i][j] = max({match, gap1, gap2});  
    }  
  }  
  
  // Backtrack to count matches
  int i = m, j = n;  
  int matches = 0;
  int totalLength = 0;
  
  while (i > 0 || j > 0) {  
    ++totalLength;
    if (i > 0 && j > 0 && dpMatrix[i][j] == dpMatrix[i - 1][j - 1] + 
    (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchPenalty)) {  
      if (sequence1[i - 1] == sequence2[j - 1]) {
        ++matches;
      }
      --i;  
      --j;  
    } else if (i > 0 && dpMatrix[i][j] == dpMatrix[i - 1][j] + gapPenalty) {  
      --i;  
    } else {  
      --j;  
    }  
  }  
  
  return static_cast<double>(matches) / totalLength;
}

// [[Rcpp::export]]
NumericMatrix calculateSimilarityMatrix(CharacterVector sequences) {
  size_t n = sequences.length();
  NumericMatrix similarityMatrix(n, n);
  
  // Calculate pairwise similarities
  for(size_t i = 0; i < n; ++i) {
    for(size_t j = i; j < n; ++j) {
      if(i == j) {
        similarityMatrix(i, j) = 1.0;  // Diagonal elements
        continue;
      }
      
      // Convert CharacterVector elements to std::string
      string seq1 = as<string>(sequences[i]);
      string seq2 = as<string>(sequences[j]);
      
      // Get similarity score directly
      double similarity = needleman_wunsch_score(seq1, seq2);
      
      // Fill both upper and lower triangles
      similarityMatrix(i, j) = similarity;
      similarityMatrix(j, i) = similarity;
    }
  }
  
  // Create numeric labels 1,2,3...
  CharacterVector labels(n);
  for(size_t i = 0; i < n; ++i) {
    labels[i] = std::to_string(i + 1);
  }
  
  // Add dimension names using numeric labels
  similarityMatrix.attr("dimnames") = List::create(labels, labels);
  
  return similarityMatrix;
}