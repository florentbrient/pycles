/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"
#include "FortranMatrix.h"

extern "C" {
  int __matrices_MOD_one_norm(int* m, int* n, int* A);
  int __matrices_MOD_one_norm_2(int* m, int* A);
  void __matrices_MOD_multiply(int* m, int* n, int* A, const int b[], int Ab[]);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
  static int one_norm(const FortranMatrix<int> *A);

  static int one_norm_2(const FortranMatrix<int> *A);

/*! \param[in] b ARRAY
 *
 *  \param[out] Ab ARRAY
*/
  static void multiply(const FortranMatrix<int> *A, const std::vector<int>* b, std::vector<int>* Ab);

};

#endif /* FORTFUNCS_H_ */
