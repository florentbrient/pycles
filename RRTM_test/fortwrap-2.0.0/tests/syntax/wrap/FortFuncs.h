/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  signed char __syntax_MOD_add_ints1(signed char* a, signed char* b);
  short __syntax_MOD_add_ints2(short* a, short* b);
  int __syntax_MOD_add_ints4(int* a, int* b);
  int __syntax_MOD_contains_arg_clash(int* a, int* b);
  int __syntax_MOD_argument_case_sensitivity(int* X);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
  static signed char add_ints1(signed char a, signed char b);

  static short add_ints2(short a, short b);

  static int add_ints4(int a, int b);

  static int contains_arg_clash(int a, int b);

  static int argument_case_sensitivity(int X);

};

#endif /* FORTFUNCS_H_ */
