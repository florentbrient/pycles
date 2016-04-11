/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  float __reals_MOD_add_floats(float* a, float* b, float* c);
  double __reals_MOD_add_doubles(double* a, double* b);
  float __reals_MOD_add_floats_lower(float* a, float* b);
  double __reals_MOD_add_doubles_lower(double* a, double* b);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
  static float add_floats(float a, float b, float c);

  static double add_doubles(double a, double b);

  static float add_floats_lower(float a, float b);

  static double add_doubles_lower(double a, double b);

};

#endif /* FORTFUNCS_H_ */
