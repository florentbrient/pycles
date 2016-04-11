/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "Object.h"

// Constructor:
Object::Object(int x) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_object(&data_ptr); // Allocate Fortran derived type
  __source_MOD_object_constructor_i(data_ptr, &x); // Fortran Constructor
  initialized = true;
}

// Constructor:
Object::Object(float f) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_object(&data_ptr); // Allocate Fortran derived type
  __source_MOD_object_constructor_f(data_ptr, &f); // Fortran Constructor
  initialized = true;
}

// Destructor:
Object::~Object() {
  __cppwrappers_MOD_deallocate_object(data_ptr); // Deallocate Fortran derived type
}

int Object::add(int i) {
  return __source_MOD_object_add(data_ptr, &i);
}

