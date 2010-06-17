#ifndef PARAMETERS_H
#define PARAMETERS_H

/*
  All run-time settings (functional weights and parameters)
  are given their own name in the enum xc_parameters.
  We can also associate descriptions and functional objects
  to these parameters.
 */

enum xc_parameters
  {
#define PARAM(name) name
#include "list_of_parameters.h"
#undef PARAM
  };

extern "C"
const char *xc_name(int param);
extern "C"
const char *xc_short_description(int param);
extern "C"
const char *xc_long_description(int param);
// Short description is a string without newlines!
void xcint_set_short_description(int param, const char *text);
// Long description should end with a final newline, and may have many lines.
void xcint_set_long_description(int param, const char *text);

class functional;
void xcint_set_functional(int param, functional *f);
functional *xcint_functional(int param);

double xcint_default(int param);

#endif
