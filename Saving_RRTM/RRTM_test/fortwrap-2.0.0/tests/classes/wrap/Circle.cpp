/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "Circle.h"

// Constructor:
Circle::Circle(int radius) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_circle(&data_ptr); // Allocate Fortran derived type
  class_data.vptr = &__classes_MOD___vtab_classes_Circle; // Get pointer to vtab
  class_data.data = data_ptr;
  __classes_MOD_circle_ctor_sub(data_ptr, &radius); // Fortran Constructor
  initialized = true;
}

// Destructor:
Circle::~Circle() {
  __cppwrappers_MOD_deallocate_circle(data_ptr); // Deallocate Fortran derived type
}

int Circle::get_area(void) {
  return __classes_MOD_circle_area(&class_data);
}

int Circle::get_diameter(void) {
  return __classes_MOD_circle_diameter(&class_data);
}

