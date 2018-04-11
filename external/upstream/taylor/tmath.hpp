#ifndef TMATH_H
#define TMATH_H

/*
  Taylor expansions of the intrinsic functions
  Ulf Ekstrom December 2010. Currently supported:
  1/x (inv)
  exp
  log
  pow
  sqrt
  cbrt (cube root)
  atan
  gauss (exp(-x^2))
  erf
  sin
  cos
  asin
  acos
  asinh
*/

#include <cassert>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _MSC_VER
#include "micromath.hpp"
#endif

// Taylor math, template style
// N is always the order of the polynomial

template <class T, int N> struct tfuns {
  static void mul(T * z, const T * x, const T * y) {
    for (int i = 0; i <= N; i++) {
      z[i] = x[0] * y[i];
      for (int j = 1; j <= i; j++)
        z[i] += x[j] * y[i - j];
    }
  }
  // z *= x
  static void multo(T * z, const T * x) {
    for (int i = N; i >= 0; i--) {
      z[i] = x[0] * z[i];
      for (int j = 1; j <= i; j++)
        z[i] += x[j] * z[i - j];
    }
  }
  // Integrates termwise, leaves x[0] undefined!
  static void integrate(T * x) {
    for (int i = N; i >= 1; i--)
      x[i] = x[i - 1] / i;
  }
  static void differentiate(T * x) {
    for (int i = 1; i <= N; i++)
      x[i - 1] = i * x[i];
    x[N] = 0;
  }
  // Set z = z + d
  static void shift(T * x, T d) {
    T dn[N + 1];
    dn[0] = 1;
    for (int i = 1; i <= N; i++)
      dn[i] = d * dn[i - 1];
    // Contributions to x[n]
    for (int n = 0; n < N; n++) {
      int fac = n + 1;
      for (int m = n + 1; m < N; m++) {
        x[n] += fac * dn[m - n] * x[m];
        fac *= m + 1;
        fac /= m - n + 1;
      }
      x[n] += fac * dn[N - n] * x[N];
    }
  }
  // assuming x[0] = 0, put sum_i f[i]x^i in f
  static void compose(T * f, const T * x) {
    switch (N) {
      case 6:
        f[6] = f[1] * x[6] + 2 * f[2] * x[5] * x[1] + 2 * f[2] * x[4] * x[2] +
               f[2] * x[3] * x[3] + 3 * f[3] * x[4] * x[1] * x[1] +
               6 * f[3] * x[3] * x[2] * x[1] + f[3] * x[2] * x[2] * x[2] +
               4 * f[4] * x[3] * x[1] * x[1] * x[1] +
               6 * f[4] * x[2] * x[2] * x[1] * x[1] +
               5 * f[5] * x[2] * x[1] * x[1] * x[1] * x[1] +
               f[6] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1];
      case 5:
        f[5] = f[1] * x[5] +
               x[1] * (2 * f[2] * x[4] +
                       x[1] * (3 * f[3] * x[3] +
                               x[1] * (4 * f[4] * x[2] + f[5] * x[1] * x[1])) +
                       3 * f[3] * x[2] * x[2]) +
               2 * f[2] * x[2] * x[3];
      case 4:
        f[4] = f[1] * x[4] +
               x[1] * (2 * f[2] * x[3] +
                       x[1] * (3 * f[3] * x[2] + f[4] * x[1] * x[1])) +
               f[2] * x[2] * x[2];
      case 3:
        f[3] = f[1] * x[3] + x[1] * (2 * f[2] * x[2] + f[3] * x[1] * x[1]);
      case 2:
        f[2] = f[1] * x[2] + f[2] * x[1] * x[1];
      case 1:
        f[1] = f[1] * x[1];
      case 0:
        break;
      default:
        assert(0 && "Unsupported order in compose()");
    }
  }
  static void stretch(T * t, T a) {
    T an = a;
    for (int i = 1; i <= N; i++) {
      t[i] *= an;
      an *= a;
    }
  }
};

// Taylor series of 1/(a+x)
template <class T, int N> static void inv_expand(T * t, const T & a) {
  assert(a != 0 && "1/(a+x) not analytic at a = 0");
  t[0] = 1 / a;
  for (int i = 1; i <= N; i++)
    t[i] = -t[i - 1] * t[0];
}

// Evaluate the taylor series of exp(x0+x)=exp(x0)*exp(x)
template <class T, int Ndeg> static void exp_expand(T * t, const T & x0) {
  T ifac = 1;
  t[0] = exp(x0);
  for (int i = 1; i <= Ndeg; i++) {
    ifac *= i;
    t[i] = t[0] / ifac;
  }
}

// Log series log(a+x) = log(1+x/a) + log(a)
template <class T, int N> static void log_expand(T * t, const T & x0) {
  assert(x0 > 0 && "log(x) not real analytic at x <= 0");
  t[0] = log(x0);
  T x0inv = 1 / x0;
  T xn = x0inv;
  for (int i = 1; i <= N; i++) {
    t[i] = (xn / double(i)) * (2 * (i & 1) - 1);
    xn *= x0inv;
  }
}

/* Use that (x0+x)^a=x0^a*(1+x/x0)^a */
template <class T, int N> static void pow_expand(T * t, T x0, T a) {
  if (x0 <= 0)
    assert(x0 > 0 && "pow(x,a) not real analytic at x <= 0");
  t[0] = pow(x0, a);
  T x0inv = 1 / x0;
  for (int i = 1; i <= N; i++)
    t[i] = t[i - 1] * x0inv * (a - i + 1) / i;
}

/* Use that (x0+x)^a=x0^a*(1+x/x0)^a */
template <class T, int N> static void sqrt_expand(T * t, const T & x0) {
  assert(x0 > 0 && "sqrt(x) not real analytic at x <= 0");
  t[0] = sqrt(x0);
  T x0inv = 1 / x0;
  for (int i = 1; i <= N; i++)
    t[i] = t[i - 1] * ((3 * x0inv) / (2 * i) - x0inv);
}

template <class T, int N> static void cbrt_expand(T * t, const T & x0) {
  assert(x0 > 0 && "pow(x,a) not real analytic at x <= 0");
  t[0] = cbrt(x0);
  T x0inv = 1 / x0;
  for (int i = 1; i <= N; i++)
    t[i] = t[i - 1] * ((4 * x0inv) / (3 * i) - x0inv);
}

// Use that d/dx atan(x) = 1/(1 + x^2),
// Taylor expand in x^2 and integrate.
template <class T, int Ndeg> static void atan_expand(T * t, T a) {
  // Calculate taylor expansion of 1/(1+a^2+x)
  T x[Ndeg + 1];
  inv_expand<T, Ndeg>(t, 1 + a * a);
  // insert x = 2*a*x + x^2
  x[0] = 0;
  if (Ndeg > 0)
    x[1] = 2 * a;
  if (Ndeg > 1)
    x[2] = 1;
  for (int i = 3; i <= Ndeg; i++)
    x[i] = 0;
  tfuns<T, Ndeg>::compose(t, x);
  // Integrate each term and set the constant
  tfuns<T, Ndeg>::integrate(t);
  t[0] = atan(a);
}

/*
   Taylor expansion of exp(-(a+x)^2) =
   exp(-a^2-2a*x)*exp(-x^2)
   Just doing a composition is unstable near 0.
 */
template <class T, int Ndeg> static void gauss_expand(T * t, const T & a) {
  exp_expand<T, Ndeg>(t, -a * a);
  tfuns<T, Ndeg>::stretch(t, -2 * a);
  T g[Ndeg + 1];
  g[0] = 1;
  for (int i = 1; i <= Ndeg; i += 2)
    g[i] = 0;
  for (int i = 1; i <= Ndeg / 2; i++)
    g[2 * i] = -g[2 * (i - 1)] / i;
  tfuns<T, Ndeg>::multo(t, g);
}

// Use that d/dx erf(x) = 2/sqrt(pi)*exp(-x^2),
// Taylor expand in x^2 and integrate.
template <class T, int Ndeg> static void erf_expand(T * t, const T & a) {
  gauss_expand<T, Ndeg>(t, a);
  for (int i = 0; i <= Ndeg; i++)
    t[i] *= 2 / sqrt(M_PI);
  tfuns<T, Ndeg>::integrate(t);
  t[0] = erf(a);
}

template <class T, int Ndeg> static void sin_expand(T * t, const T & a) {
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

template <class T, int Ndeg> static void cos_expand(T * t, const T & a) {
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

// hyperbolic arcsin function. d/dx asinh(x) = 1/sqrt(1+x^2)
// 1 + (a+x)^2 = 1+a^2 + 2ax + x^2
template <class T, int Ndeg> static void asinh_expand(T * t, const T & a) {
  T tmp[Ndeg + 1];
  tmp[0] = 1 + a * a;
  if (Ndeg > 0)
    tmp[1] = 2 * a;
  if (Ndeg > 1)
    tmp[2] = 1;
  for (int i = 3; i <= Ndeg; i++)
    tmp[i] = 0;
  pow_expand<T, Ndeg>(t, tmp[0], -0.5);
  tfuns<T, Ndeg>::compose(t, tmp);
  tfuns<T, Ndeg>::integrate(t);
  t[0] = asinh(a);
}

// arcsin function. d/dx asin(x) = 1/sqrt(1-x^2)
// 1 - (a+x)^2 = 1-a^2 - 2ax - x^2
template <class T, int Ndeg> static void asin_expand(T * t, const T & a) {
  T tmp[Ndeg + 1];
  tmp[0] = 1 - a * a;
  if (Ndeg > 0)
    tmp[1] = -2 * a;
  if (Ndeg > 1)
    tmp[2] = -1;
  for (int i = 3; i <= Ndeg; i++)
    tmp[i] = 0;
  pow_expand<T, Ndeg>(t, tmp[0], -0.5);
  tfuns<T, Ndeg>::compose(t, tmp);
  tfuns<T, Ndeg>::integrate(t);
  t[0] = asinh(a);
}

template <class T, int Ndeg> static void acosh_expand(T * t, const T & a) {
  asinh_expand(t, a);
  for (int i = 0; i < Ndeg + 1; i++)
    t[i] *= -1;
}

template <class T, int Ndeg> static void acos_expand(T * t, const T & a) {
  T tmp[Ndeg + 1];
  tmp[0] = 1 - a * a;
  if (Ndeg > 0)
    tmp[1] = -2 * a;
  if (Ndeg > 1)
    tmp[2] = -1;
  for (int i = 3; i <= Ndeg; i++)
    tmp[i] = 0;
  pow_expand<T, Ndeg>(t, tmp[0], -0.5);
  tfuns<T, Ndeg>::compose(t, tmp);
  tfuns<T, Ndeg>::integrate(t);
  for (int i = 1; i < Ndeg; i++)
    t[i] *= -1;
  t[0] = asinh(a);
}

#endif
