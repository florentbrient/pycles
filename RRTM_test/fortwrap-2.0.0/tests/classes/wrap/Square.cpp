/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "Square.h"

// Constructor:
Square::Square(int side_length) {
  data_ptr = NULL;
  __cppwrappers_MOD_allocate_square(&data_ptr); // Allocate Fortran derived type
  class_data.vptr = &__classes_MOD___vtab_classes_Square; // Get pointer to vtab
  class_data.data = data_ptr;
  __classes_MOD_square_ctor_sub(data_ptr, &side_length); // Fortran Constructor
  initialized = true;
}

// Destructor:
Square::~Square() {
  __cppwrappers_MOD_deallocate_square(data_ptr); // Deallocate Fortran derived type
}

int Square::get_area(void) {
  return __classes_MOD_square_area(&class_data);
}

