/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

int FortFuncs::array_in(const std::vector<int>* x) {
  int n = x->size();
  return __arrays_MOD_array_in(&n, &(*x)[0]);
}

void FortFuncs::array_out(std::vector<int>* x, std::vector<int>* y) {
  int nx = x->size();
  int ny = y->size();
  __arrays_MOD_array_out(&nx, &ny, &(*x)[0], &(*y)[0]);
}

int FortFuncs::inner_prod(const std::vector<int>* a, const std::vector<int>* b) {
  int n = a->size();
  return __arrays_MOD_inner_prod(&n, &(*a)[0], &(*b)[0]);
}

int FortFuncs::inner_prod_2(int n, const std::vector<int>* a, const std::vector<int>* b) {
  return __arrays_MOD_inner_prod_2(&n, &(*a)[0], &(*b)[0]);
}

