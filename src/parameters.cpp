#include <cassert>
#include "parameters.h"
#include <cstdio>
#include <cstdlib>

// TODO: Remove this and use the commone one
static void xc_die(const char *message, int code)
{
  fprintf(stderr,"XCFun fatal error %i: ",code);
  fprintf(stderr,"%s",message);
  fprintf(stderr,"\n");
  exit(-1);
}

static const char *param_symbols[XC_NR_PARAMS+1] =
  {
#define PARAM(name) #name
#include "list_of_parameters.h"
#undef PARAM
  };

static const char *param_short[XC_NR_PARAMS] = {0};
static const char *param_long[XC_NR_PARAMS] = {0};
static double param_default[XC_NR_PARAMS] = {0};

const char *xc_param_get_symbol(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);    
  return param_symbols[param];
}
const char *xc_param_get_short_description(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);
  return param_short[param];
}
const char *xc_param_get_long_description(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);
  return param_long[param];
}

double xc_param_get_default(enum xc_parameters param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);
  return param_default[param];
}

// Short description is a string without newlines!
void xcint_param_set_short_description(enum xc_parameters param, const char *text)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);
  param_short[param] = text;
}
// Long description should end with a final newline, and may have many lines.
void xcint_param_set_long_description(enum xc_parameters param, const char *text)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);
  param_long[param] = text;
}

void xcint_param_set_default(enum xc_parameters param, double value)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xc_die("Invalid parameter nr (version mismatch?)",param);
  param_default[param] = value;
}
