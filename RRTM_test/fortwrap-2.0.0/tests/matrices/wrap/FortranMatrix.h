/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTRANMATRIX_H_
#define FORTRANMATRIX_H_

#include <cstdlib>
#include <cassert>

/*! \brief A templated class for working with Matrices that store data
 *  internally in Fortran order.
 *  
 *  From C++, the data are accessed in the natural order, using
 *  <tt>x(row,column)</tt> notation, starting with base index 0
*/
template <class T>
class FortranMatrix{

  int nrows, ncols;

public:

    T *data;

/*! \brief Create a matrix with m rows and n columns
*/
  FortranMatrix(int m, int n) {
    data = NULL;
    assert(m>0 && n>0);
    nrows=m; ncols=n;
    data = (T*) calloc( m*n, sizeof(T) );
  }

  ~FortranMatrix() { if(data) free(data); }

/*! \brief Provides element access via the parentheses operator.
 *  
 *  The base index is 0
*/
  T& operator()(int i, int j) {
    assert( i>=0 && i<nrows && j>=0 && j<ncols );
    // i--; j--; // Adjust base
    return data[j*nrows+i];
  }

/*! \brief Get number of rows
*/
  inline int num_rows(void) const { return nrows; }

/*! \brief Get number of columns
*/
  inline int num_cols(void) const { return ncols; }

};

#endif /* FORTRANMATRIX_H_ */


// Local Variables:
// mode: c++
// End:
