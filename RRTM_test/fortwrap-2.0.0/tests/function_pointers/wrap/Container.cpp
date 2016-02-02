/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "Container.h"

// Constructor:
Container::Container(int (*f)(const int* a, const int* b) , int a, int b) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_container(&data_ptr); // Allocate Fortran derived type
  generic_fpointer c_pointer;
  c_pointer = (generic_fpointer) f;
  long FORT_f;
  if (f) __cppwrappers_MOD_convert_c_funcpointer(c_pointer, &FORT_f);
  __func_pointers_MOD_container_ctor(data_ptr, f ? &FORT_f : NULL, &a, &b); // Fortran Constructor
  initialized = true;
}

// Destructor:
Container::~Container() {
  __cppwrappers_MOD_deallocate_container(data_ptr); // Deallocate Fortran derived type
}

int Container::container_callf(void) {
  return __func_pointers_MOD_container_callf(data_ptr);
}

