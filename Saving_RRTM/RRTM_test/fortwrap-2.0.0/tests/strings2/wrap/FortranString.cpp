/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include "FortranString.h"

FortranString::FortranString() : length_(0), data_(NULL) {}

FortranString::FortranString(size_t length) : length_(length), data_(NULL) {
  if (length>0) data_ = (char*) calloc(length+1, sizeof(char));
}

FortranString::~FortranString() { if(data_) free(data_); }

size_t FortranString::length(void) { return length_; }

void FortranString::resize(size_t length) {
  if (data_) free(data_);
  data_ = (char*) calloc(length+1, sizeof(char));
  length_ = length;
}

void FortranString::assign(const char* s) {
  length_ = strlen(s);
  resize(length_);
  strncpy(data_, s, length_);
}

int FortranString::compare(const char* s) const {
  return strncmp(data_, s, length_);
}

char* FortranString::data(void) { return data_; }

char* FortranString::c_str(void) { return data_; }
