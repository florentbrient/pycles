/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include "InterfaceDefs.h"

extern "C" {
  int __multidim_arrays_MOD_inner_prod(int* n, const int a[], const int b[]);
  int __multidim_arrays_MOD_inner_prod_2(int* n, const int a[], const int b[]);
  void __multidim_arrays_MOD_mat_vec_mult(int* m, int* n, int* A, const int b[], int Ab[]);
  void __multidim_arrays_MOD_three_d_array_test(int* A);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
/*! \param[in] a ARRAY
 *
 *  \param[in] b ARRAY
*/
  static int inner_prod(int n, const int a[], const int b[]);

/*! \param[in] a ARRAY
 *
 *  \param[in] b ARRAY
*/
  static int inner_prod_2(int n, const int a[], const int b[]);

/*! \param[in] b ARRAY
 *
 *  \param[out] Ab ARRAY
*/
  static void mat_vec_mult(int m, int n, int* A, const int b[], int Ab[]);

  static void three_d_array_test(int* A);

};

#endif /* FORTFUNCS_H_ */
