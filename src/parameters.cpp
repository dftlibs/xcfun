#include <cassert>
#include "parameters.h"
#include "xcfun_internal.h"

static const char *param_symbols[XC_NR_PARAMS] =
  {
#define PARAM(name) #name
#include "list_of_parameters.h"
#undef PARAM
  };

static const char *param_short[XC_NR_PARAMS] = {0};
static const char *param_long[XC_NR_PARAMS] = {0};

const char *param_get_symbol(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?): ",param);    
  return param_symbols[XC_NR_PARAMS];
}
const char *param_get_short_description(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?): ",param);
  return param_short[XC_NR_PARAMS];
}
const char *param_get_long_description(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?): ",param);
  return param_long[XC_NR_PARAMS];
}

// Short description is a string without newlines!
void param_set_short_description(enum xc_parameters param, const char *text)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?): ",param);
  param_short[param] = text;
}
// Long description should end with a final newline, and may have many lines.
void param_set_long_description(enum xc_parameters param, const char *text)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?): ",param);
  param_long[param] = text;
}
