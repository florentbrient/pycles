/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include "InterfaceDefs.h"

extern "C" {
  int __c_arrays_MOD_array_in(int* n, const int x[]);
  int __c_arrays_MOD_array_in_2(const int x[]);
  void __c_arrays_MOD_array_out(int* nx, int* ny, int x[], int y[]);
  int __c_arrays_MOD_inner_prod(int* n, const int a[], const int b[]);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
/*! \param[in] x ARRAY
*/
  static int array_in(int n, const int x[]);

/*! \param[in] x ARRAY
*/
  static int array_in_2(const int x[]);

/*! \param[out] x ARRAY
 *
 *  \param[out] y ARRAY
*/
  static void array_out(int nx, int ny, int x[], int y[]);

/*! \param[in] a ARRAY
 *
 *  \param[in] b ARRAY
*/
  static int inner_prod(int n, const int a[], const int b[]);

};

#endif /* FORTFUNCS_H_ */
