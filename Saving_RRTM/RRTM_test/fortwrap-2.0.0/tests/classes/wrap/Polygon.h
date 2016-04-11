/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef POLYGON_H_
#define POLYGON_H_


#ifndef SWIG // Protect declarations from SWIG
#include <cstdlib>
#include <vector>
#include "InterfaceDefs.h"
#include "Shape.h"

// Declare external vtab data:
extern int __classes_MOD___vtab_classes_Polygon; // int is dummy data type

extern "C" {
  int __classes_MOD_polygon_num_sides(ADDRESS s);
}
#endif // SWIG

class Polygon : public Shape {

protected:
  // Polygon can not be instantiated
  Polygon() {}

public:
  virtual ~Polygon() {}

  int num_sides(void);

};

#endif /* POLYGON_H_ */
