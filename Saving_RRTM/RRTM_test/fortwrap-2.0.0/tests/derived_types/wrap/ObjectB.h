/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef OBJECTB_H_
#define OBJECTB_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"
#include "ObjectA.h"

extern "C" {
  void __cppwrappers_MOD_allocate_objectb(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_objectb(ADDRESS caddr);
  void __derivedtypes_MOD_b_ctor(ADDRESS b, ADDRESS a);
  int __derivedtypes_MOD_getax(ADDRESS b);
}
#endif // SWIG

class ObjectB {

public:
  ObjectB(ObjectA* a);
  ~ObjectB();

  int getax(void);

  ADDRESS data_ptr;

protected:
  bool initialized;
};

#endif /* OBJECTB_H_ */
