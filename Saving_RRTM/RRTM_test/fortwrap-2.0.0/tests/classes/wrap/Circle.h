/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef CIRCLE_H_
#define CIRCLE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"
#include "Shape.h"

// Declare external vtab data:
extern int __classes_MOD___vtab_classes_Circle; // int is dummy data type

extern "C" {
  void __cppwrappers_MOD_allocate_circle(ADDRESS *caddr);
  void __cppwrappers_MOD_deallocate_circle(ADDRESS caddr);
  void __classes_MOD_circle_ctor_sub(ADDRESS s, int* radius);
  int __classes_MOD_circle_area(ADDRESS s);
  int __classes_MOD_circle_diameter(ADDRESS s);
}
#endif // SWIG

/*! \brief A round object
*/
class Circle : public Shape {

public:
/*! \brief This is a user-written wrapper around the constructor above that
 *  returns a derivd type.  Currently, FortWrap does not wrap procedures
 *  that return a derived type, so if such a procedure is in the interface,
 *  a subroutine-style version should be added before running FortWrap
*/
  Circle(int radius);
  ~Circle();

/*! \brief Compute area of a circle
*/
  int get_area(void);

/*! \brief Compute area of a circle
*/
  int get_diameter(void);

};

#endif /* CIRCLE_H_ */
