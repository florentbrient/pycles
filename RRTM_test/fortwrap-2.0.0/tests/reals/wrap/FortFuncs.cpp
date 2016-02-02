/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortFuncs.h"

float FortFuncs::add_floats(float a, float b, float c) {
  return __reals_MOD_add_floats(&a, &b, &c);
}

double FortFuncs::add_doubles(double a, double b) {
  return __reals_MOD_add_doubles(&a, &b);
}

float FortFuncs::add_floats_lower(float a, float b) {
  return __reals_MOD_add_floats_lower(&a, &b);
}

double FortFuncs::add_doubles_lower(double a, double b) {
  return __reals_MOD_add_doubles_lower(&a, &b);
}

