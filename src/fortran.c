#include "xcfun.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/*
  Functions used by the Fortran interface. Do _not_ put underscores in the names!
 */
#ifndef XCFun_FORTRAN_INT
typedef int fortran_int_t;
#else
typedef XCFun_FORTRAN_INT fortran_int_t;
#endif

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

void xcint_die(const char * message, int code);

double FSYM(xcfuve) FCSYM(XCFUVE)(void) { return xcfun_version(); }

int FSYM(xcnewf) FCSYM(XCNEWF)(fortran_int_t * api_version) {
  int i;
  for (i = 0; i < MAX_FORTRAN_FUNCTIONALS; i++) {
    if (fortran_functionals[i] == 0) {
      fortran_functionals[i] = xc_new_functional_not_macro(*api_version);
      return i;
    }
  }
  return -1;
}

void FSYM(xcfree) FCSYM(XCFREE)(fortran_int_t * fun) {
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  xc_free_functional(fortran_functionals[*fun]);
  fortran_functionals[*fun] = 0;
}

void FSYM(xceval) FCSYM(XCEVAL)(fortran_int_t * fun,
                                fortran_int_t * nr_points,
                                double * first_density,
                                double * second_density,
                                double * first_result,
                                double * second_result) {
  int dpitch, rpitch;
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  /* Figure out the memory distance between points
     from the first and second pointers. */
  dpitch = second_density - first_density;
  rpitch = second_result - first_result;
  xc_eval_vec(fortran_functionals[*fun],
              *nr_points,
              first_density,
              dpitch,
              first_result,
              rpitch);
}

fortran_int_t FSYM(xcoule)
    FCSYM(XCOULE)(fortran_int_t * fun, fortran_int_t * order) {
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_output_length(fortran_functionals[*fun]);
}

#if 0
int FSYM(xcdind) FCSYM(XCDIND)(fortran_int_t *fun, const fortran_int_t *derivative)
{
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_derivative_index(fortran_functionals[*fun],derivative);
}
#endif

fortran_int_t FSYM(xcevse) FCSYM(XCTRYV)(fortran_int_t * fun,
                                         const fortran_int_t * vars,
                                         const fortran_int_t * mode,
                                         const fortran_int_t * order) {
  assert(*fun >= 0 && *fun < MAX_FORTRAN_FUNCTIONALS);
  return xc_eval_setup(fortran_functionals[*fun], *vars, *mode, *order);
}

static void str2ints(fortran_int_t ints[], fortran_int_t len, const char * s) {
  fortran_int_t i = 0;
  while (*s && i < len)
    ints[i++] = *s++;
  if (*s)
    ints[len - 1] = 0;
  else
    ints[i] = 0;
}

void FSYM(xcspla) FCSYM(XCSPLA)(fortran_int_t * text, fortran_int_t * len) {
  str2ints(text, *len, xcfun_splash());
}

int FSYM(xcsets) FCSYM(XCSETS)(fortran_int_t * fun,
                               double * value,
                               fortran_int_t * namelen,
                               const char * name) {
  int i;
  char buf[257];
  if (*namelen > 256)
    xcint_die("In xcsets_(): name string too long", *namelen);
  for (i = 0; i < *namelen; i++)
    buf[i] = name[i];
  buf[*namelen] = 0;
  return xc_set(fortran_functionals[*fun], buf, *value);
}

int FSYM(xcgets) FCSYM(XCGETS)(fortran_int_t * fun,
                               double * value,
                               fortran_int_t * namelen,
                               const char * name) {
  int i;
  char buf[257];
  if (*namelen > 256)
    xcint_die("In xcgets_(): name string too long", *namelen);
  for (i = 0; i < *namelen; i++)
    buf[i] = name[i];
  buf[*namelen] = 0;
  return xc_get(fortran_functionals[*fun], buf, value);
}

int FSYM(xcseri)
    FCSYM(XCSERI)(fortran_int_t * fun, fortran_int_t * buflen, double * buf) {
  return xc_serialize(fortran_functionals[*fun], *buflen, buf);
}
