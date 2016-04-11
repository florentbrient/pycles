/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  signed char __integers_MOD_add_ints1(signed char* a, signed char* b);
  short __integers_MOD_add_ints2(short* a, short* b);
  int __integers_MOD_add_ints4(int* a, int* b);
  long long __integers_MOD_add_ints8(long long* a, long long* b);
  signed char __integers_MOD_add_ints1_lower(signed char* a, signed char* b);
  short __integers_MOD_add_ints2_lower(short* a, short* b);
  int __integers_MOD_add_ints4_lower(int* a, int* b);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
  static signed char add_ints1(signed char a, signed char b);

  static short add_ints2(short a, short b);

  static int add_ints4(int a, int b);

  static long long add_ints8(long long a, long long b);

  static signed char add_ints1_lower(signed char a, signed char b);

  static short add_ints2_lower(short a, short b);

  static int add_ints4_lower(int a, int b);

};

#endif /* FORTFUNCS_H_ */
