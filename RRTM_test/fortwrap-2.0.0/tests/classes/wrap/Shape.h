/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef SHAPE_H_
#define SHAPE_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"

// Declare external vtab data:
extern int __classes_MOD___vtab_classes_Shape; // int is dummy data type

extern "C" {
}
#endif // SWIG

/*! \brief Base shape type
*/
class Shape {

protected:
  // Shape can not be instantiated
  Shape() {}

public:
  virtual ~Shape() {}

  virtual int get_area(void) = 0;

  ADDRESS data_ptr;
  FClassContainer class_data;

protected:
  bool initialized;
};

#endif /* SHAPE_H_ */
