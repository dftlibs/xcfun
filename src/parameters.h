#ifndef PARAMETERS_H
#define PARAMETERS_H

enum xc_parameters
  {
#define PARAM(name) name
#include "list_of_parameters.h"
#undef PARAM
  };

const char *xc_param_get_symbol(enum xc_parameters param);
const char *xc_param_get_short_description(enum xc_parameters param);
const char *xc_param_get_long_description(enum xc_parameters param);
// Short description is a string without newlines!
void xcint_param_set_short_description(enum xc_parameters param, const char *text);
// Long description should end with a final newline, and may have many lines.
void xcint_param_set_long_description(enum xc_parameters param, const char *text);

#endif
