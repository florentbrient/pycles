/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "ObjectA.h"

// Constructor:
ObjectA::ObjectA(int x) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_objecta(&data_ptr); // Allocate Fortran derived type
  __derivedtypes_MOD_a_ctor(data_ptr, &x); // Fortran Constructor
  initialized = true;
}

// Destructor:
ObjectA::~ObjectA() {
  __cppwrappers_MOD_deallocate_objecta(data_ptr); // Deallocate Fortran derived type
}

int ObjectA::getx(void) {
  return __derivedtypes_MOD_getx(data_ptr);
}

