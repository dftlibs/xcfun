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
#include "taylor.h"

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


// Rational function P(x)/Q(x), for P and Q polynomials
// of degree np and nq with coefficients pc and qc.
template<class T, class T2>
static T ratfun(const T &x, int np, const T2 pc[], int nq, const T2 qc[])
{
  return poly(x,np,pc)/poly(x,nq,qc);
}

template<class T, int Nvar, int Ndeg, class T2>
static taylor<T,Nvar,Ndeg> ratfun(const taylor<T,Nvar,Ndeg> &x, 
			   int np, const T2 pc[], int nq, const T2 qc[])
{
  taylor<T,1,Ndeg> t(x[0],0);
  return x.compose(poly(t,np,pc)/poly(t,nq,qc));
}

template<class T, int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> preexpand(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &),
				     const taylor<T,Nvar,Ndeg> &t)
{
  return t.compose(f(taylor<T,1,Ndeg>(t[0],0)));
}

template<class T, int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> preexpand(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &,
						    const T *),
			      const taylor<T,Nvar,Ndeg> &t,
			      const T *param)
{
  return t.compose(f(taylor<T,1,Ndeg>(t[0],0),param));
}

// Function of two independent variables. This is not very optimized,
// seems like a small but significant improvement for Nvar = 5, Ndeg = 2.
template<class T, int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> preexpand(taylor<T,2,Ndeg> (*f)(const taylor<T,2,Ndeg> &,
						    const taylor<T,2,Ndeg> &),
			      const taylor<T,Nvar,Ndeg> &x,
			      const taylor<T,Nvar,Ndeg> &y) 
{
  taylor<T,2,Ndeg> ft = f(taylor<T,2,Ndeg>(x[0],0),taylor<T,2,Ndeg>(y[0],1));
  taylor<T,Nvar,Ndeg> top[2];
  top[0] = x;
  top[0][0] = 0;
  top[1] = y;
  top[1][0] = 0;
  return ft.eval(top);
}


template<class T>
static T ufunc(const T &x, T &a)
{
  return pow(1+x,a)+pow(1-x,a);
}

template<class T, int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> ufunc(const taylor<T,Nvar,Ndeg> &x, const ireal_t &a)
{
  taylor<T,1,Ndeg> tmp1,tmp2;
  taylor<T,Nvar,Ndeg> res;
  pow_taylor(tmp1,1+x[0],a);
  pow_taylor(tmp2,1-x[0],a);
  for (int i=1;i<=Ndeg;i+=2)
    tmp2[i] *= -1;
  tmp1 += tmp2;
  x.compose(res,tmp1);
  return res;
}

#ifdef XCFUN_CONTRACTIONS
template<class T, int Ndeg>
static ctaylor<T,Ndeg> ufunc(const ctaylor<T,Ndeg> &x, const ireal_t &a)
{
  taylor<T,1,Ndeg> tmp1,tmp2;
  pow_taylor(tmp1,1+x[0],a);
  pow_taylor(tmp2,1-x[0],a);
  for (int i=1;i<=Ndeg;i+=2)
    tmp2[i] *= -1;
  tmp1 += tmp2;
  ctaylor<T,Ndeg> res;
  x.compose(res,tmp1.c);
  return res;
}
#endif

#endif
