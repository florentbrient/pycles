/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "ObjectB.h"

// Constructor:
ObjectB::ObjectB(ObjectA* a) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_objectb(&data_ptr); // Allocate Fortran derived type
  __derivedtypes_MOD_b_ctor(data_ptr, a->data_ptr); // Fortran Constructor
  initialized = true;
}

// Destructor:
ObjectB::~ObjectB() {
  __cppwrappers_MOD_deallocate_objectb(data_ptr); // Deallocate Fortran derived type
}

int ObjectB::getax(void) {
  return __derivedtypes_MOD_getax(data_ptr);
}

