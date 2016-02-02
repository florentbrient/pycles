/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

int FortFuncs::one_norm(const FortranMatrix<int> *A) {
  int m = A->num_rows();
  int n = A->num_cols();
  return __matrices_MOD_one_norm(&m, &n, A->data);
}

int FortFuncs::one_norm_2(const FortranMatrix<int> *A) {
  int m = A->num_rows();
  return __matrices_MOD_one_norm_2(&m, A->data);
}

void FortFuncs::multiply(const FortranMatrix<int> *A, const std::vector<int>* b, std::vector<int>* Ab) {
  int m = Ab->size();
  int n = b->size();
  __matrices_MOD_multiply(&m, &n, A->data, &(*b)[0], &(*Ab)[0]);
}

