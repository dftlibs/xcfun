#include <cassert>
#include "parameters.h"
#include "xcfun_internal.h"
#include <cstdio>
#include <cstdlib>
#if 0
// TODO: Remove this and use the commone one
static void xcint_die(const char *message, int code)
{
  fprintf(stderr,"XCFun fatal error %i: ",code);
  fprintf(stderr,"%s",message);
  fprintf(stderr,"\n");
  exit(-1);
}
#endif

static const char *param_symbols[XC_NR_PARAMS+1] =
  {
#define PARAM(name) #name
#include "list_of_parameters.h"
#undef PARAM
  };

static const char *param_short[XC_NR_PARAMS] = {0};
static const char *param_long[XC_NR_PARAMS] = {0};
static double param_default[XC_NR_PARAMS] = {0};
/* Not all of the elements will be actual functionals,
   these should be kept as null pointers. */
static functional *param_functionals[XC_NR_PARAMS] = {0};

extern "C"
const char *xc_name(int param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);    
  return param_symbols[param];
}
extern "C"
const char *xc_short_description(int param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  return param_short[param];
}
extern "C"
const char *xc_long_description(int param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  return param_long[param];
}

double xcint_default(int param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  return param_default[param];
}


// Short description is a string without newlines!
void xcint_set_short_description(int param, const char *text)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  param_short[param] = text;
}
// Long description should end with a final newline, and may have many lines.
void xcint_set_long_description(int param, const char *text)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  param_long[param] = text;
}

void xcint_set_default(int param, double value)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  param_default[param] = value;
}

// Short description is a string without newlines!
void xcint_set_functional(int param, functional *f)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  param_functionals[param] = f;
}
// Long description should end with a final newline, and may have many lines.
functional *xcint_functional(int param)
{
  if (param < 0 or param >= XC_NR_PARAMS)
      xcint_die("Invalid parameter nr (version mismatch?)",param);
  return param_functionals[param];
}
