/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef OBJECT_H_
#define OBJECT_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  void __cppwrappers_MOD_allocate_object(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_object(ADDRESS caddr);
  void __source_MOD_object_constructor_i(ADDRESS o, int* x);
  void __source_MOD_object_constructor_f(ADDRESS o, float* f);
  int __source_MOD_object_add(ADDRESS o, int* i);
}
#endif // SWIG

class Object {

public:
  Object(int x);
  Object(float f);
  ~Object();

  int add(int i);

  ADDRESS data_ptr;

protected:
  bool initialized;
};

#endif /* OBJECT_H_ */
