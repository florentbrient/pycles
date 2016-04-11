/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

int FortFuncs::callf(int (*f)(const int* a, const int* b) , int a, int b) {
  generic_fpointer c_pointer;
  c_pointer = (generic_fpointer) f;
  long FORT_f;
  if (f) __cppwrappers_MOD_convert_c_funcpointer(c_pointer, &FORT_f);
  return __func_pointers_MOD_callf(f ? &FORT_f : NULL, &a, &b);
}

