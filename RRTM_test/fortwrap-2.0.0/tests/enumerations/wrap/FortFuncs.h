/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  int __enumerations_MOD_is_red(int* x);
  int __enumerations_MOD_is_red_a(int* x);
  int __enumerations_MOD_is_green(int* x);
  int __enumerations_MOD_is_green_a(int* x);
  int __enumerations_MOD_is_blue(int* x);
  int __enumerations_MOD_is_blue_a(int* x);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
  static int is_red(int x);

  static int is_red_a(int x);

  static int is_green(int x);

  static int is_green_a(int x);

  static int is_blue(int x);

  static int is_blue_a(int x);

};

class FortConstants {
public:
  enum { RED, GREEN, BLUE };
  enum { RED_A, GREEN_A, BLUE_A };
};

#endif /* FORTFUNCS_H_ */
