#ifndef SPECMATH_H
#define SPECMATH_H

#ifdef XCFUN_NO_ERF
void xcint_die(const char *message, int code);
template<class T>
T erf(T x)
{
  xcint_die("XcFun erf called but XCFUN_NO_ERF was defined",0);
  return 0;
}
#endif

// Some math-related functions useful in many places
#include "ctaylor.h"

template<class T>
static T pow2(const T &t)
{
  return t*t;
}

template<class T>
static T pow3(const T &t)
{
  return t*t*t;
}

template<class T, class T2>
static T poly(const T &x, int ndeg, const T2 coeffs[])
{
  // Horner rule
  T res = coeffs[--ndeg];
  while (ndeg)
    {
      res *= x;
      res += coeffs[--ndeg];
    }
  return res;
}

template<class T, class S>
static T ufunc(const T &x, S a)
{
  return pow(1+x,a)+pow(1-x,a);
}

#endif
