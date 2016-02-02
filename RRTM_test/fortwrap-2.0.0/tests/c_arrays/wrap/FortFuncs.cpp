/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

int FortFuncs::array_in(int n, const int x[]) {
  return __c_arrays_MOD_array_in(&n, x);
}

int FortFuncs::array_in_2(const int x[]) {
  return __c_arrays_MOD_array_in_2(x);
}

void FortFuncs::array_out(int nx, int ny, int x[], int y[]) {
  __c_arrays_MOD_array_out(&nx, &ny, x, y);
}

int FortFuncs::inner_prod(int n, const int a[], const int b[]) {
  return __c_arrays_MOD_inner_prod(&n, a, b);
}

