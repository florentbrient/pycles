/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef SQUARE_H_
#define SQUARE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"
#include "Polygon.h"

// Declare external vtab data:
extern int __classes_MOD___vtab_classes_Square; // int is dummy data type

extern "C" {
  void __cppwrappers_MOD_allocate_square(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_square(ADDRESS caddr);
  void __classes_MOD_square_ctor_sub(ADDRESS s, int* side_length);
  int __classes_MOD_square_area(ADDRESS s);
}
#endif // SWIG

class Square : public Polygon {

public:
  Square(int side_length);
  ~Square();

  int get_area(void);

};

#endif /* SQUARE_H_ */
