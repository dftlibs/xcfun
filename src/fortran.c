#include "xcfun.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*
  Functions used by the Fortran interface. Use single underscore 
  convention. Do _not_ put underscores in the names except for
  the final one.
 */

#define FSYM(name) name##_
#define MAX_FORTRAN_FUNCTIONALS 5

static xc_functional fortran_functionals[MAX_FORTRAN_FUNCTIONALS] = {0};

double FSYM(xcfuve)(void)
{
  return xcfun_version();
}

int FSYM(xcnewf)(void)
{
  int i;
  for (i=0;i<MAX_FORTRAN_FUNCTIONALS;i++)
    {
      if (fortran_functionals[i] == 0)
	{
	  fortran_functionals[i] = xc_new_functional();
	  return i;
	}
    }
  fprintf(stderr,"Too many XC functionals (%i), check src/fortran.cpp.",MAX_FORTRAN_FUNCTIONALS);
  exit(EXIT_FAILURE);
  return -1;
}

void FSYM(xcregu)(int *fun, double *density)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_regularize_density(fortran_functionals[*fun],density);
}

void FSYM(xceval)(int *fun,  int *order, int *nr_points, double *density, double *result)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_eval(fortran_functionals[*fun], *order, *nr_points, density, result);
}

void FSYM(xcsmod)(int *fun, int *mode)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_set_mode(fortran_functionals[*fun], *mode);
}

int FSYM(xcgett)(int *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_get_type(fortran_functionals[*fun]);
}

int FSYM(xcmord)(int *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_max_order(fortran_functionals[*fun]);
}

int FSYM(xcinle)(int *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_input_length(fortran_functionals[*fun]);
}

int FSYM(xcoule)(int *fun, int *order)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_output_length(fortran_functionals[*fun],*order);
}

int FSYM(xcdind)(int *fun, const int *derivative)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_derivative_index(fortran_functionals[*fun],derivative);
}

static void str2ints(int ints[], int len, const char *s)
{
  int i = 0;
  while (*s && i < len)
    ints[i++] = *s++;
  if (*s)
    ints[len-1] = 0;
  else
    ints[i] = 0;
}

#if 0
static char *ints2str(int ints[])
{
  int len = 0;
  while (ints[len])
    len++;
  char *s = new char[len+1];
  for (int i=0;i<=len;i++)
    s[i] = ints[i];
  return s;
}
#endif

void FSYM(xcspla)(int *text, int *len)
{
  str2ints(text,*len,xcfun_splash());
}

void FSYM(xcsnam)(int *dst, int *dstlen, int *n)
{
  const char *s;
  s= xc_name(*n-1*);
  if (s)
    str2ints(dst,*dstlen,s);
  else
    dst[0] = 0;
}

void FSYM(xcssho)(int *dst, int *dstlen, int *n)
{
  const char *s;
  s = xc_short_description(*n-1);
  str2ints(dst,*dstlen,s);
}

void FSYM(xcslon)(int *dst, int *dstlen, int *n)
{
  const char *s;   
  s = xc_long_description(*n-1);
  str2ints(dst,*dstlen,s);
}

int FSYM(xcisfu)(int *n)
{
  return xc_is_functional(*n-1);
}

void FSYM(xcsets)(int *fun, int *n, double *value)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_set_param(fortran_functionals[*fun],*n-1,*value);
}

double FSYM(xcgets)(int *fun, int *n)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_get_param(fortran_functionals[*fun],*n-1);
}
