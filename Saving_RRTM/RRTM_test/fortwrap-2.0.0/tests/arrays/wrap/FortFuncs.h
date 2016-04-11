/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  int __arrays_MOD_array_in(int* n, const int x[]);
  void __arrays_MOD_array_out(int* nx, int* ny, int x[], int y[]);
  int __arrays_MOD_inner_prod(int* n, const int a[], const int b[]);
  int __arrays_MOD_inner_prod_2(int* n, const int a[], const int b[]);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
/*! \param[in] x ARRAY
*/
  static int array_in(const std::vector<int>* x);

/*! \param[out] x ARRAY
 *
 *  \param[out] y ARRAY
*/
  static void array_out(std::vector<int>* x, std::vector<int>* y);

/*! \param[in] a ARRAY
 *
 *  \param[in] b ARRAY
*/
  static int inner_prod(const std::vector<int>* a, const std::vector<int>* b);

/*! \param[in] a ARRAY
 *
 *  \param[in] b ARRAY
*/
  static int inner_prod_2(int n, const std::vector<int>* a, const std::vector<int>* b);

};

#endif /* FORTFUNCS_H_ */
