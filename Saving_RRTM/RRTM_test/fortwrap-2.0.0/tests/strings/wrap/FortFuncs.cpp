/* This source file automatically generated on 2016-01-26 using 
   FortWrap wrapper generator version 2.0.0 */

#include <cstring> // For strcpy
#include "FortFuncs.h"

void FortFuncs::string_out_literal_len(std::string *s) {
  // Declare memory to store output character data
  char s_c__[20+1];
  s_c__[20] = '\0';
  __strings_MOD_string_out_literal_len(s ? s_c__ : NULL, 20);
  if (s) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s_c__[i]==' '; i--) s_c__[i] = '\0';
    s->assign(s_c__);
  }
}

void FortFuncs::string_out_literal_len2(std::string *s) {
  // Declare memory to store output character data
  char s_c__[8+1];
  s_c__[8] = '\0';
  __strings_MOD_string_out_literal_len2(s ? s_c__ : NULL, 8);
  if (s) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=8-1; s_c__[i]==' '; i--) s_c__[i] = '\0';
    s->assign(s_c__);
  }
}

void FortFuncs::string_out_param_len(std::string *s1, std::string *s2, std::string *s3, std::string *s4) {
  // Declare memory to store output character data
  char s3_c__[20+1];
  s3_c__[20] = '\0';
  // Declare memory to store output character data
  char s2_c__[20+1];
  s2_c__[20] = '\0';
  // Declare memory to store output character data
  char s1_c__[20+1];
  s1_c__[20] = '\0';
  // Declare memory to store output character data
  char s4_c__[20+1];
  s4_c__[20] = '\0';
  __strings_MOD_string_out_param_len(s1 ? s1_c__ : NULL, s2 ? s2_c__ : NULL, s3 ? s3_c__ : NULL, s4 ? s4_c__ : NULL, 20, 20, 20, 20);
  if (s3) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s3_c__[i]==' '; i--) s3_c__[i] = '\0';
    s3->assign(s3_c__);
  }
  if (s2) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s2_c__[i]==' '; i--) s2_c__[i] = '\0';
    s2->assign(s2_c__);
  }
  if (s1) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s1_c__[i]==' '; i--) s1_c__[i] = '\0';
    s1->assign(s1_c__);
  }
  if (s4) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s4_c__[i]==' '; i--) s4_c__[i] = '\0';
    s4->assign(s4_c__);
  }
}

void FortFuncs::string_out_param_len2(std::string *s1, std::string *s2, std::string *s3, std::string *s4) {
  // Declare memory to store output character data
  char s3_c__[20+1];
  s3_c__[20] = '\0';
  // Declare memory to store output character data
  char s2_c__[20+1];
  s2_c__[20] = '\0';
  // Declare memory to store output character data
  char s1_c__[20+1];
  s1_c__[20] = '\0';
  // Declare memory to store output character data
  char s4_c__[20+1];
  s4_c__[20] = '\0';
  __strings_MOD_string_out_param_len2(s1 ? s1_c__ : NULL, s2 ? s2_c__ : NULL, s3 ? s3_c__ : NULL, s4 ? s4_c__ : NULL, 20, 20, 20, 20);
  if (s3) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s3_c__[i]==' '; i--) s3_c__[i] = '\0';
    s3->assign(s3_c__);
  }
  if (s2) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s2_c__[i]==' '; i--) s2_c__[i] = '\0';
    s2->assign(s2_c__);
  }
  if (s1) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s1_c__[i]==' '; i--) s1_c__[i] = '\0';
    s1->assign(s1_c__);
  }
  if (s4) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s4_c__[i]==' '; i--) s4_c__[i] = '\0';
    s4->assign(s4_c__);
  }
}

void FortFuncs::string_out_assumed_len(std::string *s1, int x, std::string *s2) {
  int s2_len__ = 0;
  if (s2) s2_len__ = s2->length();
  // Declare memory to store output character data
  char s2_c__[s2_len__+1];
  s2_c__[s2_len__] = '\0';
  int s1_len__ = 0;
  if (s1) s1_len__ = s1->length();
  // Declare memory to store output character data
  char s1_c__[s1_len__+1];
  s1_c__[s1_len__] = '\0';
  __strings_MOD_string_out_assumed_len(s1 ? s1_c__ : NULL, &x, s2 ? s2_c__ : NULL, s1_len__, s2_len__);
  if (s2) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=s1_len__-1; s2_c__[i]==' '; i--) s2_c__[i] = '\0';
    s2->assign(s2_c__);
  }
  if (s1) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=s1_len__-1; s1_c__[i]==' '; i--) s1_c__[i] = '\0';
    s1->assign(s1_c__);
  }
}

int FortFuncs::string_in_test(const char* s) {
  // Create C array for Fortran input string data
  char s_c__[20+1];
  if (s) {
    int i;
    strncpy(s_c__, s, 20+1); s_c__[20] = 0; // strncpy protects in case s is too long
    for (i=strlen(s_c__); i<20+1; i++) s_c__[i] = ' '; // Add whitespace for Fortran
  }
  return __strings_MOD_string_in_test(s ? s_c__ : NULL, 20);
}

int FortFuncs::string_in_cutoff(const char* s) {
  // Create C array for Fortran input string data
  char s_c__[6+1];
  if (s) {
    int i;
    strncpy(s_c__, s, 6+1); s_c__[6] = 0; // strncpy protects in case s is too long
    for (i=strlen(s_c__); i<6+1; i++) s_c__[i] = ' '; // Add whitespace for Fortran
  }
  return __strings_MOD_string_in_cutoff(s ? s_c__ : NULL, 6);
}

int FortFuncs::string_in_assumed_len(const char* s1, int x, const char* s2) {
  int s2_len__ = 0;
  if (s2) s2_len__ = strlen(s2); // Protect Optional args
  int s1_len__ = 0;
  if (s1) s1_len__ = strlen(s1); // Protect Optional args
  return __strings_MOD_string_in_assumed_len(s1, &x, s2, s1_len__, s2_len__);
}

void FortFuncs::multiple_args(int dumint, const char* s1, float dumfloat, std::string *s2) {
  // Declare memory to store output character data
  char s2_c__[20+1];
  s2_c__[20] = '\0';
  // Create C array for Fortran input string data
  char s1_c__[20+1];
  if (s1) {
    int i;
    strncpy(s1_c__, s1, 20+1); s1_c__[20] = 0; // strncpy protects in case s1 is too long
    for (i=strlen(s1_c__); i<20+1; i++) s1_c__[i] = ' '; // Add whitespace for Fortran
  }
  __strings_MOD_multiple_args(&dumint, s1 ? s1_c__ : NULL, &dumfloat, s2 ? s2_c__ : NULL, 20, 20);
  if (s2) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s2_c__[i]==' '; i--) s2_c__[i] = '\0';
    s2->assign(s2_c__);
  }
}

int FortFuncs::optional_in(const char* s) {
  // Create C array for Fortran input string data
  char s_c__[20+1];
  if (s) {
    int i;
    strncpy(s_c__, s, 20+1); s_c__[20] = 0; // strncpy protects in case s is too long
    for (i=strlen(s_c__); i<20+1; i++) s_c__[i] = ' '; // Add whitespace for Fortran
  }
  return __strings_MOD_optional_in(s ? s_c__ : NULL, 20);
}

int FortFuncs::optional_in_assumed(const char* s) {
  int s_len__ = 0;
  if (s) s_len__ = strlen(s); // Protect Optional args
  return __strings_MOD_optional_in_assumed(s, s_len__);
}

int FortFuncs::optional_out(std::string *s) {
  // Declare memory to store output character data
  char s_c__[20+1];
  s_c__[20] = '\0';
  int __retval = __strings_MOD_optional_out(s ? s_c__ : NULL, 20);
  if (s) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=20-1; s_c__[i]==' '; i--) s_c__[i] = '\0';
    s->assign(s_c__);
  }
  return __retval;
}

int FortFuncs::optional_out_assumed(std::string *s) {
  int s_len__ = 0;
  if (s) s_len__ = s->length();
  // Declare memory to store output character data
  char s_c__[s_len__+1];
  s_c__[s_len__] = '\0';
  int __retval = __strings_MOD_optional_out_assumed(s ? s_c__ : NULL, s_len__);
  if (s) {
    // Trim trailing whitespace and assign character array to string:
    for (int i=s_len__-1; s_c__[i]==' '; i--) s_c__[i] = '\0';
    s->assign(s_c__);
  }
  return __retval;
}

