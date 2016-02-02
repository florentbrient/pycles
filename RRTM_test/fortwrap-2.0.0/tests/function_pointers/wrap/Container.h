/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef CONTAINER_H_
#define CONTAINER_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  void __cppwrappers_MOD_allocate_container(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_container(ADDRESS caddr);
  void __func_pointers_MOD_container_ctor(ADDRESS c, void* f, int* a, int* b);
  int __func_pointers_MOD_container_callf(ADDRESS c);
}
#endif // SWIG

class Container {

public:
  Container(int (*f)(const int* a, const int* b) , int a, int b);
  ~Container();

  int container_callf(void);

  ADDRESS data_ptr;

protected:
  bool initialized;
};

#endif /* CONTAINER_H_ */
