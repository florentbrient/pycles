/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

int FortFuncs::add_allopt(const int* a, const int* b, const int* c) {
  return __opt_args_MOD_add_allopt(a, b, c);
}

int FortFuncs::add_mixed(int a, int b, const int* c, const int* d) {
  return __opt_args_MOD_add_mixed(&a, &b, c, d);
}

