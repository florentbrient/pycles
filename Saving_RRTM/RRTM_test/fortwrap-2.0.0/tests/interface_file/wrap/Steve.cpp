/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "Steve.h"

// Constructor:
Steve::Steve(int i) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_object_to_rename(&data_ptr); // Allocate Fortran derived type
  __source_MOD_object2_constructor(data_ptr, &i); // Fortran Constructor
  initialized = true;
}

// Destructor:
Steve::~Steve() {
  __cppwrappers_MOD_deallocate_object_to_rename(data_ptr); // Deallocate Fortran derived type
}

int Steve::val(void) {
  return __source_MOD_object2_val(data_ptr);
}

