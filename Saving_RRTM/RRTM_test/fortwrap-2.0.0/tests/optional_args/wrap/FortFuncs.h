/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  int __opt_args_MOD_add_allopt(const int* a, const int* b, const int* c);
  int __opt_args_MOD_add_mixed(int* a, int* b, const int* c, const int* d);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
/*! \param[in] a OPTIONAL
 *
 *  \param[in] b OPTIONAL
 *
 *  \param[in] c OPTIONAL
*/
  static int add_allopt(const int* a=NULL, const int* b=NULL, const int* c=NULL);

/*! \param[in] c OPTIONAL
 *
 *  \param[in] d OPTIONAL
*/
  static int add_mixed(int a, int b, const int* c=NULL, const int* d=NULL);

};

#endif /* FORTFUNCS_H_ */
