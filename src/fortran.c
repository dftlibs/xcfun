#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "xcfun.h"
/*
  Functions used by the Fortran interface. Use single underscore 
  convention. Do _not_ put underscores in the names except for
  the final one.
 */

#ifdef FTN_UPPERCASE
#define FSYM(name)
#define FCSYM(name) name
#elif FTN_UPPERCASE_UNDERSCORE
#define FSYM(name)
#define FCSYM(name) name##_
#else
#define FSYM(name) name##_
#define FCSYM(name)
#endif

#define MAX_FORTRAN_FUNCTIONALS 5

static xc_functional fortran_functionals[MAX_FORTRAN_FUNCTIONALS] = {0};

double FSYM(xcfuve) FCSYM(XCFUVE)(void)
{
  return xcfun_version();
}

int FSYM(xcnewf) FCSYM(XCNEWF)(void)
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
  return -1;
}

void FSYM(xcfree)FCSYM(XCFREE)(int *fun)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_free_functional(fortran_functionals[*fun]);
  fortran_functionals[*fun] = 0;
}

void FSYM(xceval)FCSYM(XCEVAL)(int *fun, 
		  int *nr_points, 
		  double *first_density,
		  double *second_density,
		  double *first_result,
		  double *second_result)
{
  int dpitch, rpitch;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
/* Figure out the memory distance between points
   from the first and second pointers. */
  dpitch = second_density - first_density;
  rpitch = second_result - first_result;
  xc_eval_vec(fortran_functionals[*fun], *nr_points, 
	      first_density, dpitch,
	      first_result, rpitch);
}

int FSYM(xcoule)FCSYM(XCOULE)(int *fun, int *order)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_output_length(fortran_functionals[*fun]);
}

#if 0
int FSYM(xcdind)FCSYM(XCDIND)(int *fun, const int *derivative)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_derivative_index(fortran_functionals[*fun],derivative);
}
#endif

int FSYM(xcevse)FCSYM(XCTRYV)(int *fun, const int *vars, const int *mode, const int *order)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_eval_setup(fortran_functionals[*fun], *vars, *mode, *order);
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

void FSYM(xcspla)FCSYM(XCSPLA)(int *text, int *len)
{
  str2ints(text,*len,xcfun_splash());
}

#if 0
void FSYM(xcsnam)FCSYM(XCSNAM)(int *dst, int *dstlen, int *n)
{
  const char *s;
  s= xc_name(*n-1);
  if (s)
    str2ints(dst,*dstlen,s);
  else
    dst[0] = 0;
}

void FSYM(xcssho)FCSYM(XCSSHO)(int *dst, int *dstlen, int *n)
{
  const char *s;
  s = xc_short_description(*n-1);
  str2ints(dst,*dstlen,s);
}

void FSYM(xcslon)FCSYM(XCSLON)(int *dst, int *dstlen, int *n)
{
  const char *s;   
  s = xc_long_description(*n-1);
  str2ints(dst,*dstlen,s);
}

int FSYM(xcisfu)FCSYM(XCISFU)(int *n)
{
  return xc_is_functional(*n-1);
}

#endif

void FSYM(xcsets)FCSYM(XCSETS)(int *fun, int *n, double *value)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_set(fortran_functionals[*fun],*n,*value);
}

double FSYM(xcgets)FCSYM(XCGETS)(int *fun, int *n)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_get(fortran_functionals[*fun],*n);
}

