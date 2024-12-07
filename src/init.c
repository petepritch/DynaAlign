#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

// External declaration
extern SEXP minhash_similarity_matrix(SEXP, SEXP, SEXP);
extern SEXP calculateSimilarityMatrix(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef callMethods[] = {
  {"minhash_similarity_matrix", (DL_FUNC) &minhash_similarity_matrix, 3},
  {"calculateSimilarityMatrix", (DL_FUNC) &calculateSimilarityMatrix, 4},
  {NULL, NULL, 0}
};

void R_init_DynaAlign(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}