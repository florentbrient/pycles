/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

int FortFuncs::inner_prod(int n, const int a[], const int b[]) {
  return __multidim_arrays_MOD_inner_prod(&n, a, b);
}

int FortFuncs::inner_prod_2(int n, const int a[], const int b[]) {
  return __multidim_arrays_MOD_inner_prod_2(&n, a, b);
}

void FortFuncs::mat_vec_mult(int m, int n, int* A, const int b[], int Ab[]) {
  __multidim_arrays_MOD_mat_vec_mult(&m, &n, A, b, Ab);
}

void FortFuncs::three_d_array_test(int* A) {
  __multidim_arrays_MOD_three_d_array_test(A);
}

