#include <Rcpp.h>
#include <string>
#include <vector>
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
List needleman_wunsch(const string &sequence1, const string &sequence2, int matchScore, int mismatchPenalty, int gapPenalty) {  
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
      int match = dpMatrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchPenalty);  
      int gap1 = dpMatrix[i - 1][j] + gapPenalty;  
      int gap2 = dpMatrix[i][j - 1] + gapPenalty;  
      
      dpMatrix[i][j] = max({match, gap1, gap2});  
    }  
  }  
  
  // Backtrack to find the optimal alignment  
  int i = m, j = n;  
  string align1, align2;  
  
  while (i > 0 || j > 0) {  
    if (i > 0 && j > 0 && dpMatrix[i][j] == dpMatrix[i - 1][j - 1] + (sequence1[i - 1] == sequence2[j - 1] ? matchScore : mismatchPenalty)) {  
      align1 = sequence1[i - 1] + align1;  
      align2 = sequence2[j - 1] + align2;  
      --i;  
      --j;  
    } else if (i > 0 && dpMatrix[i][j] == dpMatrix[i - 1][j] + gapPenalty) {  
      align1 = sequence1[i - 1] + align1;  
      align2 = '-' + align2;  
      --i;  
    } else {  
      align1 = '-' + align1;  
      align2 = sequence2[j - 1] + align2;  
      --j;  
    }  
  }  
  return List::create(Named("alignment1") = align1,
                      Named("alignment2") = align2);
}