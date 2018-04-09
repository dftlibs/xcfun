// By Ulf Ekstrom March-June 2009.
// Elementary functions. For inclusion in taylor.h only!

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _MSC_VER
#include "micromath.hpp"
#endif

// Taylor series of 1/(a+x)
template <class T, int N> static void inv_taylor(taylor<T, 1, N> & t, const T & a) {
  assert(a != 0 && "1/(a+x) not analytic at a = 0");
  t[0] = 1 / a;
  for (int i = 1; i <= N; i++)
    t[i] = -t[i - 1] * t[0];
}

// Reciprocal taylor series by the Newton method.
// a/t = 1 + O(x^Ndeg+1)
// t0 = 1/a[0]
// t1 = 1/a[0]*(2-a/a[0])
// tn+1 = -tn*(a*tn-2)
template <class T, int Nvar, int Ndeg>
static void taylor_reciprocal(taylor<T, Nvar, Ndeg> & t,
                              const taylor<T, Nvar, Ndeg> & a) {
  assert(a != 0 && "1/(a+x) not analytic at a = 0");
  t[0] = 1 / a[0];
  if (Ndeg == 0)
    return;
  T c = -t[0] * t[0];
  for (int i = 1; i < t.size; i++)
    t[i] = c * a[i];

  // need only log2(Ndeg) number of terms
  for (int i = 1; i <= (int)log2(Ndeg + 1); i++) {
    taylor<T, Nvar, Ndeg> tmp = a;
    tmp *= t;
    tmp[0] -= 2;
    tmp *= t;
    for (int j = 0; j < t.size; j++)
      t[j] = -tmp[j];
  }
}

template <class T, int Nvar, int Ndeg, class S>
static taylor<T, Nvar, Ndeg> operator/(const S & x,
                                       const taylor<T, Nvar, Ndeg> & t) {
#ifdef TAYLOR_LOGGING
  if (taylor_logging)
    cout << "operator/ S templated. x = " << x << endl;
#endif
  taylor<T, 1, Ndeg> tmp;
  inv_taylor(tmp, t[0]);
  tmp *= x;
  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> operator/(const taylor<T, Nvar, Ndeg> & t,
                                       const T & x) {
#ifdef TAYLOR_LOGGING
  if (taylor_logging)
    cout << "operator/ by scalar. x = " << x << endl;
#endif
  taylor<T, Nvar, Ndeg> tmp = t;
  tmp *= 1 / x;
  return tmp;
}

template <class T, int Nvar, int Ndeg, class S>
static taylor<T, Nvar, Ndeg> operator/(const taylor<T, Nvar, Ndeg> & t,
                                       const S & x) {
#ifdef TAYLOR_LOGGING
  if (taylor_logging)
    cout << "operator/ by scalar templated. x = " << x << endl;
#endif
  taylor<T, Nvar, Ndeg> tmp = t;
  tmp *= 1 / T(x);
  return tmp;
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> operator/(const taylor<T, Nvar, Ndeg> & t1,
                                       const taylor<T, Nvar, Ndeg> & t2) {
  taylor<T, 1, Ndeg> tmp;
  inv_taylor(tmp, t2[0]);
  taylor<T, Nvar, Ndeg> res;
  t2.compose(res, tmp);
  res *= t1;
  return res;
}

// Evaluate the taylor series of exp(x0+x)=exp(x0)*exp(x)
template <class T, int Ndeg>
static void exp_taylor(taylor<T, 1, Ndeg> & t, const T & x0) {
  T ifac = 1;
  t[0] = exp(x0);
  for (int i = 1; i <= Ndeg; i++) {
    ifac *= i;
    t[i] = t[0] / ifac;
  }
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> exp(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  exp_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

// exp(x)-1, but accurate for small x
template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> expm1(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  exp_taylor(tmp, t[0]);

  // Only constant value is affected by the cancellation
  if (fabs(t[0]) > 1e-3)
    tmp[0] -= 1;
  else
    tmp[0] = 2 * exp(t[0] / 2) * sinh(t[0] / 2);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

// Log series log(a+x) = log(1+x/a) + log(a)
template <class T, int N>
static void log_taylor_old(taylor<T, 1, N> & t, const T & x0) {
  //  assert(x0 != T(0) && "log(x) not analytic at x = 0");
  t[0] = log(x0);
  T xn = x0;
  for (int i = 1; i <= N; i++) {
    t[i] = (2 * (i & 1) - 1) / (i * xn);
    xn *= x0;
  }
}

// Log series log(a+x) = log(1+x/a) + log(a)
template <class T, int N> static void log_taylor(taylor<T, 1, N> & t, const T & x0) {
  assert(x0 > 0 && "log(x) not real analytic at x <= 0");
  t[0] = log(x0);
  T x0inv = 1 / x0;
  T xn = x0inv;
  for (int i = 1; i <= N; i++) {
    t[i] = (xn / double(i)) * (2 * (i & 1) - 1);
    xn *= x0inv;
  }
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> log(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  log_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

/* Use that (x0+x)^a=x0^a*(1+x/x0)^a */
template <class S, class T, int N>
static void pow_taylor(taylor<T, 1, N> & t, const T & x0, const S & a) {
  if (x0 <= 0)
    assert(x0 > 0 && "pow(x,a) not real analytic at x <= 0");
  t[0] = pow(x0, a);
  T x0inv = 1 / x0;
  for (int i = 1; i <= N; i++)
    t[i] = t[i - 1] * x0inv * (a - i + 1) / i;
}

// We need this version with double a argument to prevent truncation
// to int.
template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> pow(const taylor<T, Nvar, Ndeg> & t, const double & a) {
  taylor<T, 1, Ndeg> tmp;
  pow_taylor(tmp, t[0], a);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

template <class T, int Nvar, int Ndeg, class S>
static taylor<T, Nvar, Ndeg> pow(const taylor<T, Nvar, Ndeg> & t, const S & a) {
  taylor<T, 1, Ndeg> tmp;
  pow_taylor(tmp, t[0], T(a));

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

#if 0
// This might be slightly dangerous to enable this, since
// we can construct and a from an int.
// when the exponent is also a taylor expansion use
// t^a = exp(a*log(t))
template<class T,int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> pow(const taylor<T,Nvar,Ndeg> &t, 
			const taylor<T,Nvar,Ndeg> &a)
{
  return exp(a*log(t));
}
#endif

/* Use that (x0+x)^a=x0^a*(1+x/x0)^a */
template <class T, int N>
static void sqrt_taylor(taylor<T, 1, N> & t, const T & x0) {
  assert(x0 > 0 && "sqrt(x) not real analytic at x <= 0");
  t[0] = sqrt(x0);
  T x0inv = 1 / x0;
  for (int i = 1; i <= N; i++)
    t[i] = t[i - 1] * ((3 * x0inv) / (2 * i) - x0inv);
  // x0inv*(3.0/(2*i) - 1.0);
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> sqrt(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  sqrt_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

// This is used only when there is no native cbrt.
template <class T> inline T cbrt(T x) { return pow(x, 1.0 / 3.0); }

template <class T, int N>
static void cbrt_taylor(taylor<T, 1, N> & t, const T & x0) {
  assert(x0 > 0 && "pow(x,a) not real analytic at x <= 0");
  t[0] = cbrt(x0);
  T x0inv = 1 / x0;
  for (int i = 1; i <= N; i++)
    t[i] = t[i - 1] * ((4 * x0inv) / (3 * i) - x0inv);
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> cbrt(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  cbrt_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

// Integer exponent version is analytical at t[0] = 0
// This function gets priority over the normal pow
// when the exponent is an integer, but does not force
// conversion to integer.
template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> pow(const taylor<T, Nvar, Ndeg> & t, int n) {
  if (n > 0) {
    taylor<T, Nvar, Ndeg> res = t;
    while (n-- > 1)
      res *= t;
    return res;
  } else if (n < 0) {
    return 1 / pow(t, -n);
  } else {
    taylor<T, Nvar, Ndeg> res(1);
    return res;
  }
}

// Use that d/dx atan(x) = 1/(1 + x^2),
// Taylor expand in x^2 and integrate.
template <class T, int Ndeg>
static void atan_taylor(taylor<T, 1, Ndeg> & t, const T & a) {
  // Calculate taylor expansion of 1/(1+a^2+x)
  taylor<T, 1, Ndeg> invt, x;
  inv_taylor(invt, 1 + a * a);
  // insert x = 2*a*x + x^2
  x = 0;
  if (Ndeg > 0)
    x[1] = 2 * a;
  if (Ndeg > 1)
    x[2] = 1;
  x.compose0(t, invt);
  // Integrate each term and set the constant
  t.integrate();
  t[0] = atan(a);
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> atan(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  atan_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

/*
   Taylor expansion of exp(-(a+x)^2) =
   exp(-a^2-2a*x)*exp(-x^2)
 */
template <class T, int Ndeg>
static void gauss_taylor(taylor<T, 1, Ndeg> & t, const T & a) {
  exp_taylor(t, -a * a);
  t.stretch(-2 * a);
  taylor<T, 1, Ndeg> g;
  g = 1;
  for (int i = 1; i <= Ndeg / 2; i++)
    g[2 * i] = -g[2 * (i - 1)] / i;
  t *= g;
}
#ifdef WITH_QD
// Does not work afaik
/* // QD does not provide an error function, so here is a slow and */
/* // not so smart Taylor expansion. */
/* static inline qd_real erf(const qd_real &x) */
/* { */
/*   qd_real sum = 0; */
/*   qd_real x2 = x*x, term = x; */
/*   int k = 0; */
/*   if (fabs(x) > 12) // erf(12) = 1 - 1e-64 */
/*     return 1; */
/*   sum += term; */
/*   if (fabs(x) > 1e-64) */
/*     { */
/*       double magn = erf(to_double(x)); // order of magnitude of answer. */
/*       while (fabs(term/magn) > 1e-64) */
/* 	{ */
/* 	  k++; */
/* 	  term *= qd_real(-(2*k-1)*x2)/(k*(2*k+1)); */
/* 	  sum += term; */
/* 	} */
/*     } */
/*   return sum*2/sqrt(qd_real::_pi); */
/* } */
#endif

// Use that d/dx erf(x) = 2/sqrt(pi)*exp(-x^2),
// Taylor expand in x^2 and integrate.
template <class T, int Ndeg>
static void erf_taylor(taylor<T, 1, Ndeg> & t, const T & a) {
  gauss_taylor(t, a);
  t *= 2 / sqrt(M_PI);
  // Integrate each term and set the constant
  t.integrate();
  t[0] = erf(a);
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> erf(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  erf_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

template <class T, int Ndeg>
static void sin_taylor(taylor<T, 1, Ndeg> & t, const T & a) {
  if (Ndeg > 0) {
    T s = sin(a), c = cos(a), fac = 1;
    for (int i = 0; 2 * i < Ndeg; i++) {
      t[2 * i] = fac * s;
      fac /= (2 * i + 1);
      t[2 * i + 1] = fac * c;
      fac /= -(2 * i + 2);
    }
    if (Ndeg % 2 == 0)
      t[Ndeg] = s * fac;
  } else {
    t[0] = sin(a);
  }
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> sin(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  sin_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

template <class T, int Ndeg>
static void cos_taylor(taylor<T, 1, Ndeg> & t, const T & a) {
  if (Ndeg > 0) {
    T s = sin(a), c = cos(a), fac = 1;
    for (int i = 0; 2 * i < Ndeg; i++) {
      t[2 * i] = fac * c;
      fac /= -(2 * i + 1);
      t[2 * i + 1] = fac * s;
      fac /= (2 * i + 2);
    }
    if (Ndeg % 2 == 0)
      t[Ndeg] = c * fac;
  } else {
    t[0] = cos(a);
  }
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> cos(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  cos_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

// hyperbolic arcsin function. d/dx asinh(x) = 1/sqrt(1+x^2)
// 1 + (a+x)^2 = 1+a^2 + 2ax + x^2
template <class T, int Ndeg>
static void asinh_taylor(taylor<T, 1, Ndeg> & t, const T & a) {
  taylor<T, 1, Ndeg> tmp(1 + a * a);
  if (Ndeg > 0)
    tmp[1] = 2 * a;
  if (Ndeg > 1)
    tmp[2] = 1;
  t = pow(tmp, -0.5);
  t.integrate();
  t[0] = asinh(a);
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> asinh(const taylor<T, Nvar, Ndeg> & t) {
  taylor<T, 1, Ndeg> tmp;
  asinh_taylor(tmp, t[0]);

  taylor<T, Nvar, Ndeg> res;
  t.compose(res, tmp);
  return res;
}

template <class T, int Ndeg> static void sinc_taylor_at0(taylor<T, 1, Ndeg> & t) {
  t[0] = 1;
  T fac = 1;
  for (int i = 1; i <= Ndeg / 2; i++) {
    fac /= (2 * i);
    t[2 * i - 1] = 0;
    fac /= -(2 * i + 1);
    t[2 * i] = fac;
  }
  if (Ndeg % 2 == 1)
    t[Ndeg] = 0;
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> sinc(const taylor<T, Nvar, Ndeg> & t) {
  if (fabs(t[0]) < 1e-3) {
    taylor<T, 1, Ndeg + 8> tmp;
    sinc_taylor_at0(tmp);
    taylor<T, 1, Ndeg> tmp2;
    T dx = t[0];
    tmp.shift(tmp2, &dx);
    taylor<T, Nvar, Ndeg> res;
    t.compose(res, tmp2);
    return res;
  } else {
    return sin(t) / t;
  }
}

/*
  The original function is unstable for small t[0] values similarly to
  the Boys function. Use an [8,8] Pade approximation when |t[0]| is
  small. This works less well but still ok in single precision.
 */
template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> sqrtx_asinh_sqrtx(const taylor<T, Nvar, Ndeg> & t) {
  assert(t[0] > -0.5);
  // Coefficients of an [8,8] Pade approximation at x = 0
  static const T P[] = {0,
                        3.510921856028398e3,
                        1.23624388373212e4,
                        1.734847003883674e4,
                        1.235072285222234e4,
                        4.691117148130619e3,
                        9.119186273274577e2,
                        7.815848629220836e1,
                        1.96088643023654e0};
  static const T Q[] = {3.510921856028398e3,
                        1.29475924799926e4,
                        1.924308297963337e4,
                        1.474357149568687e4,
                        6.176496729255528e3,
                        1.379806958043824e3,
                        1.471833349002349e2,
                        5.666278232986776e0,
                        2.865104054302032e-2};
  if (fabs(t[0]) < 0.5) {
    // Shift polys
    assert(Ndeg <= 8);
    taylor<T, 1, Ndeg> pp, pq, ppade;
    reinterpret_cast<const taylor<T, 1, 8> *>(P)->shift(pp, &t[0]);
    reinterpret_cast<const taylor<T, 1, 8> *>(Q)->shift(pq, &t[0]);
    ppade = pp / pq;
    taylor<T, Nvar, Ndeg> res;
    t.compose(res, ppade);
    return res;
  } else {
    // This is the unstable form
    taylor<T, Nvar, Ndeg> s = sqrt(t);
    return s * asinh(s);
  }
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> min(const taylor<T, Nvar, Ndeg> & a,
                                 const taylor<T, Nvar, Ndeg> & b) {
  if (a[0] <= b[0])
    return a;
  else
    return b;
}

template <class T, int Nvar, int Ndeg>
static taylor<T, Nvar, Ndeg> max(const taylor<T, Nvar, Ndeg> & a,
                                 const taylor<T, Nvar, Ndeg> & b) {
  if (a[0] > b[0])
    return a;
  else
    return b;
}
