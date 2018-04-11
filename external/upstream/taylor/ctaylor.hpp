#ifndef CTAYLOR_H
#define CTAYLOR_H
#include <cassert>
#include <cmath>
#include <cstdio>

// sparse(scalar) is not cost effective at second order,
// maybe at third? Can also use NAN in second
// coefficient to make scalars.
//#define CTAYLOR_SPARSE

#define CNST 0 // avoid defining CONST
#define VAR0 1
#define VAR1 2
#define VAR2 4
#define VAR3 8
#define VAR4 16
#define VAR5 32
#define VAR6 64
#define VAR7 128

/*
  ctaylor is a class of tensored first order polynomials,
  giving a total of 2^Nvar coefficients, with the highest
  order term being of order Nvar. C is for cube. Coefficients
  are stored with the first half of the terms not depending
  on the last variable, while the second half does. This
  is repeated recursively. For example with Nvar=3:
  1 x y xy z zx zy zxy
 */

#define POW2(N) (1 << (N))

/*
  Recursion on the ctaylor data.

  P_n = P_n-1 + x_nR_n-1
  P_n*Q_n = P_n-1*Q_n-1 + P_n-1*Sn-1 + R_n-1*Q_n

 */
template <class T, int Nvar> struct ctaylor_rec {
  // Add x*y to dst
  static void mul(T * dst, const T * x, const T * y) {
    ctaylor_rec<T, Nvar - 1>::mul(dst, x, y);
    ctaylor_rec<T, Nvar - 1>::mul(dst + POW2(Nvar - 1), x + POW2(Nvar - 1), y);
    ctaylor_rec<T, Nvar - 1>::mul(dst + POW2(Nvar - 1), x, y + POW2(Nvar - 1));
  }
  // Add set dst = x*y to dst
  static void mul_set(T * dst, const T * x, const T * y) {
    ctaylor_rec<T, Nvar - 1>::mul_set(dst, x, y);
    ctaylor_rec<T, Nvar - 1>::mul_set(dst + POW2(Nvar - 1), x + POW2(Nvar - 1), y);
    ctaylor_rec<T, Nvar - 1>::mul(dst + POW2(Nvar - 1), x, y + POW2(Nvar - 1));
  }
  // dst = dst * y
  static void multo(T * dst, const T * y) {
    ctaylor_rec<T, Nvar - 1>::multo(dst + POW2(Nvar - 1), y);
    ctaylor_rec<T, Nvar - 1>::mul(dst + POW2(Nvar - 1), dst, y + POW2(Nvar - 1));
    ctaylor_rec<T, Nvar - 1>::multo(dst, y);
  }
  // dst = dst * (y - y[0])
  static void multo_skipconst(T * dst, const T * y) {
    ctaylor_rec<T, Nvar - 1>::multo_skipconst(dst + POW2(Nvar - 1), y);
    ctaylor_rec<T, Nvar - 1>::mul(dst + POW2(Nvar - 1), dst, y + POW2(Nvar - 1));
    ctaylor_rec<T, Nvar - 1>::multo_skipconst(dst, y);
  }

  // Evaluate the polynomial at x (an array of length Nvar)
  static T eval(const T * ct, const T * x) {
    return ctaylor_rec<T, Nvar - 1>::eval(ct, x) +
           x[Nvar - 1] * ctaylor_rec<T, Nvar - 1>::eval(ct + POW2(Nvar - 1), x);
  }
  /* Put sum_i coeff[i]*(x - x[0])^i in res,
     used when evaluating analytical functions of this */
  static void compose(T * res, const T * x, const T coeff[]) {
    res[0] = coeff[Nvar];
    for (int i = 1; i < POW2(Nvar); i++)
      res[i] = 0;
    for (int i = Nvar - 1; i >= 0; i--) {
      ctaylor_rec<T, Nvar>::multo_skipconst(res, x);
      res[0] += coeff[i];
    }
  }
};

template <class T> struct ctaylor_rec<T, 0> {
  static void mul(T * dst, const T * x, const T * y) { dst[0] += x[0] * y[0]; }
  static void mul_set(T * dst, const T * x, const T * y) { dst[0] = x[0] * y[0]; }
  static void multo(T * dst, const T * y) { dst[0] *= y[0]; }
  static void multo_skipconst(T * dst, const T * y) { dst[0] = 0; }
  static T eval(const T * ct, const T * x) { return ct[0]; }
  static void compose(T * res, const T * x, const T coeff[]) { res[0] = coeff[0]; }
};

template <class T> struct ctaylor_rec<T, 1> {
  static void mul(T * dst, const T * x, const T * y) {
    dst[0] += x[0] * y[0];
    dst[1] += x[0] * y[1] + x[1] * y[0];
  }
  static void mul_set(T * dst, const T * x, const T * y) {
    dst[0] = x[0] * y[0];
    dst[1] = x[0] * y[1] + x[1] * y[0];
  }
  static void multo(T * dst, const T * y) {
    dst[1] = dst[1] * y[0] + dst[0] * y[1];
    dst[0] *= y[0];
  }
  static void multo_skipconst(T * dst, const T * y) {
    dst[1] = dst[0] * y[1];
    dst[0] = 0;
  }
  static T eval(const T * ct, const T * x) { return ct[0] + x[0] * ct[1]; }
  static void compose(T * res, const T * x, const T coeff[]) {
    res[0] = coeff[0];
    res[1] = coeff[1] * x[1];
  }
};

template <class T> struct ctaylor_rec<T, 2> {
  static void mul(T * dst, const T * x, const T * y) {
    dst[0] += x[0] * y[0];
    dst[1] += x[0] * y[1] + x[1] * y[0];
    dst[2] += x[0] * y[2] + x[2] * y[0];
    dst[3] += x[0] * y[3] + x[3] * y[0] + x[1] * y[2] + x[2] * y[1];
  }
  static void mul_set(T * dst, const T * x, const T * y) {
    dst[0] = x[0] * y[0];
    dst[1] = x[0] * y[1] + x[1] * y[0];
    dst[2] = x[0] * y[2] + x[2] * y[0];
    dst[3] = x[0] * y[3] + x[3] * y[0] + x[1] * y[2] + x[2] * y[1];
  }
  static void multo(T * dst, const T * y) {
    dst[3] = dst[0] * y[3] + dst[3] * y[0] + dst[1] * y[2] + dst[2] * y[1];
    dst[2] = dst[0] * y[2] + dst[2] * y[0];
    dst[1] = dst[0] * y[1] + dst[1] * y[0];
    dst[0] = dst[0] * y[0];
  }
  static void multo_skipconst(T * dst, const T * y) {
    dst[3] = dst[0] * y[3] + dst[1] * y[2] + dst[2] * y[1];
    dst[2] = dst[0] * y[2];
    dst[1] = dst[0] * y[1];
    dst[0] = 0;
  }
  static T eval(const T * ct, const T * x) {
    return ct[0] + x[0] * ct[1] + x[1] * ct[2] + x[1] * x[0] * ct[3];
  }
  static void compose(T * res, const T * x, const T coeff[]) {
    res[0] = coeff[0];
    res[1] = coeff[1] * x[1];
    res[2] = coeff[1] * x[2];
    res[3] = coeff[1] * x[3] + 2 * x[1] * x[2] * coeff[2];
  }
};

template <class T, int Nvar> struct ctaylor {
  enum { size = POW2(Nvar) };
  T c[size];
#ifdef CTAYLOR_SPARSE
  int isscalar; // If isscalar the only c[0] is defined.
#endif
  ctaylor() {
#ifdef CTAYLOR_SPARSE
    isscalar = 0;
#endif
  }
  ctaylor(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
    isscalar = t.isscalar;
    if (t.isscalar) {
      c[0] = t.c[0];
    } else {
      for (int i = 0; i < POW2(Nvar); i++)
        c[i] = t.c[i];
    }
#else
    for (int i = 0; i < POW2(Nvar); i++)
      c[i] = t.c[i];
#endif
  }
  ctaylor(const T & c0) {
    c[0] = c0;
#ifdef CTAYLOR_SPARSE
    isscalar = 1;
#else
    for (int i = 1; i < POW2(Nvar); i++)
      c[i] = 0;
#endif
  }
  ctaylor(const T & c0, int var) {
    c[0] = c0;
    for (int i = 1; i < POW2(Nvar); i++)
      c[i] = 0;
    assert(var >= 0);
    assert(Nvar > var);
    c[var] = 1;
#ifdef CTAYLOR_SPARSE
    isscalar = 0;
#endif
  }
  ctaylor(const T & c0, int var, const T & varval) {
    c[0] = c0;
    for (int i = 1; i < POW2(Nvar); i++)
      c[i] = 0;
    assert(var >= 0);
    assert(Nvar > var);
    c[var] = varval;
#ifdef CTAYLOR_SPARSE
    isscalar = 0;
#endif
  }
  template <typename S> ctaylor<T, Nvar> & operator=(const S & c0) {
    c[0] = c0;
#ifdef CTAYLOR_SPARSE
    isscalar = 1;
#else
    for (int i = 1; i < POW2(Nvar); i++)
      c[i] = 0;
#endif
    return *this;
  }
  ctaylor<T, Nvar> & operator=(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
    isscalar = t.isscalar;
    if (t.isscalar)
      c[0] = t.c[0];
    else
      for (int i = 0; i < POW2(Nvar); i++)
        c[i] = t.c[i];
#else
    for (int i = 0; i < POW2(Nvar); i++)
      c[i] = t.c[i];
#endif
    return *this;
  }
  void make_nonscalar() {
#ifdef CTAYLOR_SPARSE
    isscalar = 0;
    for (int i = 1; i < POW2(Nvar); i++)
      c[i] = 0;
#endif
  }
  void set(int i, T value) {
    assert(i >= 0);
    assert(POW2(Nvar) > i);
#ifdef CTAYLOR_SPARSE
    if (i > 0 && isscalar)
      make_nonscalar();
#endif
    c[i] = value;
  }
  T get(int i) const {
    assert(i >= 0);
    assert(POW2(Nvar) > i);
#ifdef CTAYLOR_SPARSE
    if (i > 0 && isscalar)
      return 0;
    else
      return c[i];
#else
    return c[i];
#endif
  }

  ctaylor<T, Nvar> operator-(void) const {
    ctaylor<T, Nvar> res;
#ifdef CTAYLOR_SPARSE
    res.isscalar = isscalar;
    if (isscalar)
      res.c[0] = -c[0];
    else
      for (int i = 0; i < POW2(Nvar); i++)
        res.c[i] = -c[i];
#else
    for (int i = 0; i < POW2(Nvar); i++)
      res.c[i] = -c[i];
#endif
    return res;
  }
  void operator-=(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
    c[0] -= t.c[0];
    if (!t.isscalar) {
      if (isscalar)
        for (int i = 1; i < POW2(Nvar); i++)
          c[i] = -t.c[i];
      else
        for (int i = 1; i < POW2(Nvar); i++)
          c[i] -= t.c[i];
    }
    isscalar &= t.isscalar;
#else
    for (int i = 0; i < POW2(Nvar); i++)
      c[i] -= t.c[i];
#endif
  }
  void operator+=(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
    c[0] += t.c[0];
    if (!t.isscalar) {
      if (isscalar)
        for (int i = 1; i < POW2(Nvar); i++)
          c[i] = t.c[i];
      else
        for (int i = 1; i < POW2(Nvar); i++)
          c[i] += t.c[i];
    }
    isscalar &= t.isscalar;
#else
    for (int i = 0; i < POW2(Nvar); i++)
      c[i] += t.c[i];
#endif
  }
  template <class S> void operator-=(const S & x) { c[0] -= x; }
  template <class S> void operator+=(const S & x) { c[0] += x; }
  template <class S> void operator*=(const S & scale) {
#ifdef CTAYLOR_SPARSE
    if (isscalar)
      c[0] *= scale;
    else
      for (int i = 0; i < POW2(Nvar); i++)
        c[i] *= scale;
#else
    for (int i = 0; i < POW2(Nvar); i++)
      c[i] *= scale;
#endif
  }
  template <class S> void operator/=(const S & scale) {
#ifdef CTAYLOR_SPARSE
    if (isscalar)
      c[0] /= scale;
    else
      for (int i = 0; i < POW2(Nvar); i++)
        c[i] /= scale;
#else
    for (int i = 0; i < POW2(Nvar); i++)
      c[i] /= scale;
#endif
  }
  void operator/=(const ctaylor<T, Nvar> & t) { *this = *this / t; }
  void operator*=(const ctaylor<T, Nvar> & t) {
#ifdef CTAYLOR_SPARSE
    if (t.isscalar)
      if (isscalar)
        c[0] *= t.c[0];
      else
        for (int i = 0; i < POW2(Nvar); i++)
          c[i] *= t.c[0];
    else if (isscalar)
      for (int i = POW2(Nvar) - 1; i >= 0; i--)
        c[i] = t.c[0] * c[0];
    else
      ctaylor_rec<T, Nvar>::multo(c, t.c);
    isscalar &= t.isscalar;
#else
    ctaylor_rec<T, Nvar>::multo(c, t.c);
#endif
  }
};

template <class T, int Nvar>
static ctaylor<T, Nvar> operator*(const ctaylor<T, Nvar> & t1,
                                  const ctaylor<T, Nvar> & t2) {
  ctaylor<T, Nvar> tmp;
#ifdef CTAYLOR_SPARSE
  if (t2.isscalar)
    if (t1.isscalar)
      tmp.c[0] = t1.c[0] * t2.c[0];
    else
      for (int i = 0; i < POW2(Nvar); i++)
        tmp.c[i] = t1.c[i] * t2.c[0];
  else if (t1.isscalar)
    for (int i = 0; i < POW2(Nvar); i++)
      tmp.c[i] = t1.c[0] * t2.c[i];
  else
    ctaylor_rec<T, Nvar>::mul_set(tmp.c, t1.c, t2.c);
  tmp.isscalar = t1.isscalar & t2.isscalar;
#else
  ctaylor_rec<T, Nvar>::mul_set(tmp.c, t1.c, t2.c);
#endif
  return tmp;
}

// <> comparisons are taken to mean comparing the constant
// coefficient. This makes the transition from numbers to
// taylor objects easier.
template <class S, class T, int Nvar>
static bool operator<(const S & x, const ctaylor<T, Nvar> & t) {
  return x < t.c[0];
}

template <class S, class T, int Nvar>
static bool operator<(const ctaylor<T, Nvar> & t, const S & x) {
  return t.c[0] < x;
}

template <class T, int Nvar>
static bool operator<(const ctaylor<T, Nvar> & t1, const ctaylor<T, Nvar> & t2) {
  return t1.c[0] < t2.c[0];
}

template <class S, class T, int Nvar>
static bool operator>(const S & x, const ctaylor<T, Nvar> & t) {
  return x > t.c[0];
}

template <class S, class T, int Nvar>
static bool operator>(const ctaylor<T, Nvar> & t, const S & x) {
  return t.c[0] > x;
}

template <class T, int Nvar>
static bool operator>(const ctaylor<T, Nvar> & t1, const ctaylor<T, Nvar> & t2) {
  return t1.c[0] > t2.c[0];
}

template <class S, class T, int Nvar>
static bool operator!=(const ctaylor<T, Nvar> & t, const S & x) {
  return t.c[0] != x;
}

template <class T, int Nvar, class S>
static ctaylor<T, Nvar> operator*(const S & x, const ctaylor<T, Nvar> & t) {
  ctaylor<T, Nvar> tmp;
#ifdef CTAYLOR_SPARSE
  tmp.isscalar = t.isscalar;
  if (tmp.isscalar)
    tmp.c[0] = x * t.c[0];
  else
    for (int i = 0; i < POW2(Nvar); i++)
      tmp.c[i] = x * t.c[i];
#else
  for (int i = 0; i < POW2(Nvar); i++)
    tmp.c[i] = x * t.c[i];
#endif
  return tmp;
}

template <class T, int Nvar, class S>
static ctaylor<T, Nvar> operator*(const ctaylor<T, Nvar> & t, const S & x) {
  ctaylor<T, Nvar> tmp;
#ifdef CTAYLOR_SPARSE
  tmp.isscalar = t.isscalar;
  if (tmp.isscalar)
    tmp.c[0] = x * t.c[0];
  else
    for (int i = 0; i < POW2(Nvar); i++)
      tmp.c[i] = x * t.c[i];
#else
  for (int i = 0; i < POW2(Nvar); i++)
    tmp.c[i] = x * t.c[i];
#endif
  return tmp;
}

template <class T, int Nvar, class S>
static ctaylor<T, Nvar> operator+(const S & x, const ctaylor<T, Nvar> & t) {
  ctaylor<T, Nvar> tmp = t;
  tmp.c[0] += x;
  return tmp;
}

template <class T, int Nvar, class S>
static ctaylor<T, Nvar> operator+(const ctaylor<T, Nvar> & t, const S & x) {
  ctaylor<T, Nvar> tmp = t;
  tmp.c[0] += x;
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator+(const ctaylor<T, Nvar> & t1,
                                  const ctaylor<T, Nvar> & t2) {
  ctaylor<T, Nvar> tmp = t1;
  tmp += t2;
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator-(const T & x, const ctaylor<T, Nvar> & t) {
  ctaylor<T, Nvar> tmp = -t;
  tmp.c[0] += x;
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator-(const ctaylor<T, Nvar> & t, const T & x) {
  ctaylor<T, Nvar> tmp = t;
  tmp.c[0] -= x;
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator-(int x, const ctaylor<T, Nvar> & t) {
  ctaylor<T, Nvar> tmp = -t;
  tmp.c[0] += x;
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator-(const ctaylor<T, Nvar> & t, int x) {
  ctaylor<T, Nvar> tmp = t;
  tmp.c[0] -= x;
  return tmp;
}

template <class T, int Nvar>
static ctaylor<T, Nvar> operator-(const ctaylor<T, Nvar> & t1,
                                  const ctaylor<T, Nvar> & t2) {
  ctaylor<T, Nvar> tmp = t1;
  tmp -= t2;
  return tmp;
}

#include "ctaylor_math.hpp"

#ifdef CTAYLOR_CXXIO

#include <iostream>

template <class num, int Nvar>
static std::ostream & operator<<(std::ostream & stream,
                                 const ctaylor<num, Nvar> & t) {
  stream << "{" << t[0];
  for (int i = 1; i < POW2(Nvar); i++)
    stream << ", " << t[i];
  stream << "}";
  return stream;
}

#endif

#endif
