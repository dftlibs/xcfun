// For inclusion in ctaylor.h only!
#include "tmath.hpp"

template <class T, int Nvar, class S>
static ctaylor<T, Nvar> operator/(const S & x, const ctaylor<T, Nvar> & t) {
  ctaylor<T, Nvar> res;
#ifdef CTAYLOR_SPARSE
  res.isscalar = t.isscalar;
  if (res.isscalar) {
    res.c[0] = x / t.c[0];
  } else {
    T tmp[Nvar + 1];
    inv_expand<T, Nvar>(tmp, t.c[0]);
    for (int i = 0; i <= Nvar; i++)
      tmp[i] *= x;
    ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  }
#else
  T tmp[Nvar + 1];
  inv_expand<T, Nvar>(tmp, t.c[0]);
  for (int i = 0; i <= Nvar; i++)
    tmp[i] *= x;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
#endif
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator/(const ctaylor<T, Nvar> & t1,
                                  const ctaylor<T, Nvar> & t2) {
  ctaylor<T, Nvar> res;
#ifdef CTAYLOR_SPARSE
  if (t1.isscalar)
    return t1.c[0] / t2;
  else if (t2.isscalar)
    return t1 / t2.c[0];
  res.isscalar = 0;
#endif
  T tmp[Nvar + 1];
  inv_expand<T, Nvar>(tmp, t2.c[0]);
  ctaylor_rec<T, Nvar>::compose(res.c, t2.c, tmp);
  res *= t1;
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator/(const ctaylor<T, Nvar> & t, const T & x) {
  ctaylor<T, Nvar> tmp = t;
  tmp *= 1 / x;
  return tmp;
}

template <class T, int Nvar, class S>
static ctaylor<T, Nvar> operator/(const ctaylor<T, Nvar> & t, const S & x) {
  ctaylor<T, Nvar> tmp = t;
  tmp *= 1 / T(x);
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> abs(const ctaylor<T, Nvar> & t) {
  if (t.c[0] < 0)
    return -t;
  else
    return t;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> exp(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(exp(t.c[0]));
#endif
  T tmp[Nvar + 1];
  exp_expand<T, Nvar>(tmp, t.c[0]);
  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

// exp(x)-1, but accurate for small x
template <class T, int Nvar>
static ctaylor<T, Nvar> expm1(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(2 * exp(t.c[0] / 2) * sinh(t.c[0] / 2));
#endif
  T tmp[Nvar + 1];
  exp_expand<T, Nvar>(tmp, t.c[0]);

  // Only constant value is affected by the cancellation
  if (fabs(t.c[0]) > 1e-3)
    tmp[0] -= 1;
  else
    tmp[0] = 2 * exp(t.c[0] / 2) * sinh(t.c[0] / 2);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> log(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(log(t.c[0]));
#endif
  T tmp[Nvar + 1];
  log_expand<T, Nvar>(tmp, t.c[0]);
  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

// We need this version with double a argument to prevent truncation
// to int.
template <class T, int Nvar>
static ctaylor<T, Nvar> pow(const ctaylor<T, Nvar> & t, const double & a) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(pow(t.c[0], a));
#endif
  T tmp[Nvar + 1];
  pow_expand<T, Nvar>(tmp, t.c[0], a);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> sqrt(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(sqrt(t.c[0]));
#endif
  T tmp[Nvar + 1];
  sqrt_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> cbrt(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(cbrt(t.c[0]));
#endif
  T tmp[Nvar + 1];
  cbrt_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

// Integer exponent version is analytical at t[0] = 0
// This function gets priority over the normal pow
// when the exponent is an integer, but does not force
// conversion to integer.
template <class T, int Nvar>
static ctaylor<T, Nvar> pow(const ctaylor<T, Nvar> & t, int n) {
  if (n > 0) {
    ctaylor<T, Nvar> res = t;
    while (n-- > 1)
      res *= t;
    return res;
  } else if (n < 0) {
    return pow(t, double(n));
  } else {
    ctaylor<T, Nvar> res(1);
    return res;
  }
}

template <class T, int Nvar>
static ctaylor<T, Nvar> atan(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(atan(t.c[0]));
#endif
  T tmp[Nvar + 1];
  atan_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> erf(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(erf(t.c[0]));
#endif
  T tmp[Nvar + 1];
  erf_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> sin(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(sin(t.c[0]));
#endif
  T tmp[Nvar + 1];
  sin_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> cos(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(cos(t.c[0]));
#endif
  T tmp[Nvar + 1];
  cos_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> asin(const ctaylor<T, Nvar> & t) {
  T tmp[Nvar + 1];
  asin_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> acos(const ctaylor<T, Nvar> & t) {
  T tmp[Nvar + 1];
  acos_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> asinh(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
  if (t.isscalar)
    return ctaylor<T, Nvar>(asinh(t.c[0]));
#endif
  T tmp[Nvar + 1];
  asinh_expand<T, Nvar>(tmp, t.c[0]);

  ctaylor<T, Nvar> res;
  ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
  return res;
}

/*
  The original function is unstable for small t[0] values similarly to
  the Boys function. Use an [8,8] Pade approximation when |t[0]| is
  small. This works less well but still ok in single precision.
 */
template <class T, int Nvar>
static ctaylor<T, Nvar> sqrtx_asinh_sqrtx(const ctaylor<T, Nvar> & t) {
  assert(t.c[0] > -0.5);
#ifdef CTAYLOR_SPARSE
  if (t.isscalar) {
    T sqrtx = sqrt(t.c[0]);
    return ctaylor<T, Nvar>(sqrtx * asinh(sqrtx));
  }
#endif
// Coefficients of an [8,8] Pade approximation at x = 0
#define ASINH_TABSIZE 9
  static const T P[ASINH_TABSIZE] = {0,
                                     3.510921856028398e3,
                                     1.23624388373212e4,
                                     1.734847003883674e4,
                                     1.235072285222234e4,
                                     4.691117148130619e3,
                                     9.119186273274577e2,
                                     7.815848629220836e1,
                                     1.96088643023654e0};
  static const T Q[ASINH_TABSIZE] = {3.510921856028398e3,
                                     1.29475924799926e4,
                                     1.924308297963337e4,
                                     1.474357149568687e4,
                                     6.176496729255528e3,
                                     1.379806958043824e3,
                                     1.471833349002349e2,
                                     5.666278232986776e0,
                                     2.865104054302032e-2};
  if (fabs(t.c[0]) < 0.5) {
    // Shift polys, divide and compose
    assert(Nvar < ASINH_TABSIZE);
    T tmp[Nvar + 1], pq[9];
    for (int i = 0; i < ASINH_TABSIZE; i++)
      pq[i] = Q[i];
    tfuns<T, ASINH_TABSIZE - 1>::shift(pq, t.c[0]);
    inv_expand<T, Nvar>(tmp, pq[0]);
    tfuns<T, Nvar>::compose(tmp, pq);
    for (int i = 0; i < ASINH_TABSIZE; i++)
      pq[i] = P[i];
    tfuns<T, ASINH_TABSIZE - 1>::shift(pq, t.c[0]);
    tfuns<T, Nvar>::multo(tmp, pq);
    ctaylor<T, Nvar> res;
    ctaylor_rec<T, Nvar>::compose(res.c, t.c, tmp);
    return res;
  } else {
    // This is the unstable form
    ctaylor<T, Nvar> s = sqrt(t);
    return s * asinh(s);
  }
}

template <class T, int Nvar>
static ctaylor<T, Nvar> min(const ctaylor<T, Nvar> & a, const ctaylor<T, Nvar> & b) {
  if (a <= b)
    return a;
  else
    return b;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> max(const ctaylor<T, Nvar> & a, const ctaylor<T, Nvar> & b) {
  if (a > b)
    return a;
  else
    return b;
}
