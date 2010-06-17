#include "xcfun_internal.h"

#define MAX_FORTRAN_FUNCTIONALS 5

static xc_functional fortran_functionals[MAX_FORTRAN_FUNCTIONALS] = {0};

extern "C"
double xcfuve_(void)
{
  return xcfun_version();
}

extern "C"
int xcnewf_(void)
{
  for (int i=0;i<MAX_FORTRAN_FUNCTIONALS;i++)
    {
      if (fortran_functionals[i] == 0)
	{
	  fortran_functionals[i] = xc_new_functional();
	  return i;
	}
    }
  xc_die("Too many XC functionals, check src/fortran.cpp.",MAX_FORTRAN_FUNCTIONALS);
  return -1;
}

extern "C"
void xcregu_(int *fun, double *density)
{
  fortran_functionals[*fun]->regularize_density(density);
}

extern "C"
void xceval_(int *fun, double *result, int *order, double *density)
{
  fortran_functionals[*fun]->eval(result,*order,density);
}

extern "C"
void xcsmod_(int *fun, int *mode)
{
  fortran_functionals[*fun]->set_mode(*mode);
}

extern "C"
int xcgett_(int *fun)
{
  return fortran_functionals[*fun]->get_type();
}

extern "C"
int xcmord_(int *fun)
{
  return fortran_functionals[*fun]->get_max_order();
}

extern "C"
int xcinle_(int *fun)
{
  return fortran_functionals[*fun]->input_length();
}

extern "C"
int xcoule_(int *fun, int *order)
{
  return fortran_functionals[*fun]->output_length(*order);
}

extern "C"
int xcdind_(int *fun, const int *derivative)
{
  return fortran_functionals[*fun]->derivative_index(derivative);
}

static void str2ints(int ints[], int len, const char *s)
{
  int i = 0;
  while (*s and i < len)
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

extern "C"
void xcspla_(int *text, int *len)
{
  str2ints(text,*len,xcfun_splash());
}

extern "C"
void xcsnam_(int *fun, int *dst, int *dstlen, int *n)
{
  const char *s = xc_name(*n);
  if (s)
    str2ints(dst,*dstlen,s);
  else
    dst[0] = 0;
}

extern "C"
void xcssho_(int *fun, int *dst, int *dstlen, int *n)
{
  const char *s = xc_short_description(*n);
  str2ints(dst,*dstlen,s);
}

extern "C"
void xcslon_(int *fun, int *dst, int *dstlen, int *n)
{
  const char *s = xc_long_description(*n);
  str2ints(dst,*dstlen,s);
}

extern "C"
int xcisfu_(int *fun, int *n)
{
  return xc_is_functional(*n);
}

extern "C"
void xcsets_(int *fun, int *n, double *value)
{
  xc_set(fortran_functionals[*fun],*n,*value);
}

extern "C"
double xcgets_(int *fun, int *n)
{
  return xc_get(fortran_functionals[*fun],*n);
}


