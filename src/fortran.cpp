#include "xcfun_internal.h"

static array<xc_functional *> fortran_functionals;

extern "C"
double xcfuve_(void)
{
  return xcfun_version();
}




extern "C"
int xcnewf_(void)
{
  fortran_functionals.push_back(new xc_functional);
  return fortran_functionals.size()-1;
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

extern "C"
void xcspla_(int *text, int *len)
{
  str2ints(text,*len,xcfun_splash());
}

extern "C"
void xcsnam_(int *fun, int *dst, int *dstlen, int *n)
{
  const char *s = fortran_functionals[*fun]->setting_name(*n);
  if (s)
    str2ints(dst,*dstlen,s);
  else
    dst[0] = 0;
}

extern "C"
void xcssho_(int *fun, int *dst, int *dstlen, int *name)
{
  char *n = ints2str(name);
  const char *s = fortran_functionals[*fun]->setting_short_description(n);
  str2ints(dst,*dstlen,s);
  delete[] n;
}

extern "C"
void xcslon_(int *fun, int *dst, int *dstlen, int *name)
{
  char *n = ints2str(name);
  const char *s = fortran_functionals[*fun]->setting_long_description(n);
  str2ints(dst,*dstlen,s);
  delete[] n;
}

extern "C"
int xcisfu_(int *fun, int *name)
{
  char *n = ints2str(name);
  int res = fortran_functionals[*fun]->is_functional(n);
  delete[] n;
  return res;
}

extern "C"
int xcsets_(int *fun, int *name, double *value)
{
  char *n = ints2str(name);
  int res = 
    fortran_functionals[*fun]->set_setting(n,*value);
  delete[] n;
  return res;
}

extern "C"
double xcgets_(int *fun, int *name)
{
  char *n = ints2str(name);
  double d = fortran_functionals[*fun]->get_setting(n);
  delete[] n;
  return d;
}

extern "C"
int xcisse_(int *fun, int *name)
{
  char *n = ints2str(name);
  int res = fortran_functionals[*fun]->is_set(n);
  delete[] n;
  return res;
}


