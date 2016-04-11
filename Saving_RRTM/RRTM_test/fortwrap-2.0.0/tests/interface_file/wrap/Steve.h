/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef STEVE_H_
#define STEVE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

extern "C" {
  void __cppwrappers_MOD_allocate_object_to_rename(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_object_to_rename(ADDRESS caddr);
  void __source_MOD_object2_constructor(ADDRESS o, int* i);
  int __source_MOD_object2_val(ADDRESS o);
}
#endif // SWIG

class Steve {

public:
  Steve(int i);
  ~Steve();

  int val(void);

  ADDRESS data_ptr;

protected:
  bool initialized;
};

#endif /* STEVE_H_ */
