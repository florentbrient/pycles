/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTFUNCS_H_
#define FORTFUNCS_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"
#include "FortranString.h"

extern "C" {
  void __strings2_MOD_string_out_literal_len(char* s, int s_len__);
  void __strings2_MOD_string_out_literal_len2(char* s, int s_len__);
  void __strings2_MOD_string_out_param_len(char* s1, char* s2, char* s3, char* s4, int s1_len__, int s2_len__, int s3_len__, int s4_len__);
  void __strings2_MOD_string_out_param_len2(char* s1, char* s2, char* s3, char* s4, int s1_len__, int s2_len__, int s3_len__, int s4_len__);
  void __strings2_MOD_string_out_assumed_len(char* s1, int* x, char* s2, int s1_len__, int s2_len__);
  int __strings2_MOD_string_in_test(const char* s, int s_len__);
  int __strings2_MOD_string_in_cutoff(const char* s, int s_len__);
  int __strings2_MOD_string_in_assumed_len(const char* s1, int* x, const char* s2, int s1_len__, int s2_len__);
  void __strings2_MOD_multiple_args(int* dumint, const char* s1, float* dumfloat, char* s2, int s1_len__, int s2_len__);
  int __strings2_MOD_optional_in(const char* s, int s_len__);
  int __strings2_MOD_optional_in_assumed(const char* s, int s_len__);
  int __strings2_MOD_optional_out(char* s, int s_len__);
  int __strings2_MOD_optional_out_assumed(char* s, int s_len__);
}
#endif // SWIG

/*! \brief Wrapper class for Fortran routines that do not operate on a derived type
*/
class FortFuncs {

public:
  static void string_out_literal_len(FortranString *s);

  static void string_out_literal_len2(FortranString *s);

  static void string_out_param_len(FortranString *s1, FortranString *s2, FortranString *s3, FortranString *s4);

  static void string_out_param_len2(FortranString *s1, FortranString *s2, FortranString *s3, FortranString *s4);

  static void string_out_assumed_len(FortranString *s1, int x, FortranString *s2);

  static int string_in_test(const char* s);

/*! \brief Make sure conversion gets handled properly when dummy length is less
 *  than actual length
*/
  static int string_in_cutoff(const char* s);

  static int string_in_assumed_len(const char* s1, int x, const char* s2);

  static void multiple_args(int dumint, const char* s1, float dumfloat, FortranString *s2);

/*! \param[in] s OPTIONAL
*/
  static int optional_in(const char* s=NULL);

/*! \param[in] s OPTIONAL
*/
  static int optional_in_assumed(const char* s=NULL);

/*! \param[out] s OPTIONAL
*/
  static int optional_out(FortranString *s=NULL);

/*! \param[out] s OPTIONAL
*/
  static int optional_out_assumed(FortranString *s=NULL);

};

#endif /* FORTFUNCS_H_ */
