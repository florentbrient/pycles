/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef OBJECTA_H_
#define OBJECTA_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  void __cppwrappers_MOD_allocate_objecta(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_objecta(ADDRESS caddr);
  void __derivedtypes_MOD_a_ctor(ADDRESS a, int* x);
  int __derivedtypes_MOD_getx(ADDRESS a);
}
#endif // SWIG

class ObjectA {

public:
  ObjectA(int x);
  ~ObjectA();

  int getx(void);

  ADDRESS data_ptr;

protected:
  bool initialized;
};

#endif /* OBJECTA_H_ */
