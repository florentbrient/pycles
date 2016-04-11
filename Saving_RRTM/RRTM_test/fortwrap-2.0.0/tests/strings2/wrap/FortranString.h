/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#ifndef FORTRANSTRING_H_
#define FORTRANSTRING_H_

#include <cstdlib>
#include <cstring>

/*! \brief Simple wrapper class for handling dynamic allocation of string data
 *  
 *  Used by FortWrap to handle wrapping of character output arguments.  Emulates some of the basic functionality of std::string, but can be used as a way to remove dependency on std::string and avoid C++ library conflicts in some cases
*/
class FortranString{
  size_t length_;
  char* data_;
  
 public:

  FortranString();

  FortranString(size_t length);

  ~FortranString();

  size_t length(void);

  void resize(size_t length);

  void assign(const char* s);

  int compare(const char* s) const;

  char* data(void);

  char* c_str(void);
};


#endif /* FORTRANSTRING_H_ */


// Local Variables:
// mode: c++
// End:
