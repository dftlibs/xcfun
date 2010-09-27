#ifndef CTAYLOR_H
#define CTAYLOR_H
#include <cmath>
#include <cassert>

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
template<class T, int Nvar>
struct ctaylor_rec
{
  // Add x*y to dst
  static void mul(T *dst, const T *x, const T *y)
  {
    ctaylor_rec<T,Nvar-1>::mul(dst,x,y);
    ctaylor_rec<T,Nvar-1>::mul(dst+POW2(Nvar-1),x+POW2(Nvar-1),y);
    ctaylor_rec<T,Nvar-1>::mul(dst+POW2(Nvar-1),x,y+POW2(Nvar-1));
  }
  // dst = dst * y
  static void multo(T *dst, const T *y)
  {
    ctaylor_rec<T,Nvar-1>::multo(dst+POW2(Nvar-1),y);
    ctaylor_rec<T,Nvar-1>::mul(dst+POW2(Nvar-1),dst,y+POW2(Nvar-1));
    ctaylor_rec<T,Nvar-1>::multo(dst,y);
  }
  // dst = dst * (y - y[0])
  static void multo_skipconst(T *dst, const T *y)
  {
    ctaylor_rec<T,Nvar-1>::multo_skipconst(dst+POW2(Nvar-1),y);
    ctaylor_rec<T,Nvar-1>::mul(dst+POW2(Nvar-1),dst,y+POW2(Nvar-1));
    ctaylor_rec<T,Nvar-1>::multo_skipconst(dst,y);
  }

  // Evaluate the polynomial at x (an array of length Nvar)
  static T eval(const T *ct, const T *x)
  {
    return ctaylor_rec<T,Nvar-1>::eval(ct,x) 
      + x[Nvar-1]*ctaylor_rec<T,Nvar-1>::eval(ct+POW2(Nvar-1),x); 
  }
};

template<class T>
struct ctaylor_rec<T,0>
{
  static void mul(T *dst, const T *x, const T *y)
  {
    *dst += *x * *y;
  }
  static void multo(T *dst, const T *y)
  {
    *dst *= *y;
  }
  static void multo_skipconst(T *dst, const T *y)
  {
    *dst = 0;
  }
  static T eval(const T *ct, const T *x)
  {
    return ct[0];
  }
};

template<class T, int Nvar>
struct ctaylor
{
  T c[POW2(Nvar)];
  ctaylor(void) {}
  ctaylor(const T &c0)
  {
    c[0] = c0;
    for (int i=1;i<POW2(Nvar);i++)
      c[i] = 0;
  }
  ctaylor(const T &c0, int var)
  {
    c[0] = c0;
    for (int i=1;i<POW2(Nvar);i++)
      c[i] = 0;
    assert(var>=0);
    assert(Nvar>var);
    c[POW2(var)] = 1;
  }
  ctaylor(const T &c0, int var, const T &varval)
  {
    c[0] = c0;
    for (int i=1;i<POW2(Nvar);i++)
      c[i] = 0;
    assert(var>=0);
    assert(Nvar>var);
    c[POW2(var)] = varval;
  }
  template<typename S>
  ctaylor<T,Nvar> &operator=(const S& c0)
  {
    c[0] = c0;
    for (int i=1;i<POW2(Nvar);i++)
      c[i] = 0;
    return *this;
  }
  const T &operator[](int i) const
  {
    assert(i>=0);
    assert(POW2(Nvar)>i);
    return c[i];
  }
  T &operator[](int i) 
  {
    assert(i>=0);
    assert(POW2(Nvar)>i);
    return c[i];
  }
  ctaylor<T,Nvar> operator-(void) const
  {
    ctaylor<T,Nvar> res = *this;
    for (int i=0;i<POW2(Nvar);i++)
      res[i] = -res[i];
    return res;
  }
  void operator-=(const ctaylor<T,Nvar>& t)
  {
    for (int i=0;i<POW2(Nvar);i++)
      c[i] -= t.c[i];
  }
  void operator+=(const ctaylor<T,Nvar>& t)
  {
    for (int i=0;i<POW2(Nvar);i++)
      c[i] += t.c[i];
  }
  template<class S>
  void operator-=(const S& x)
  {
    c[0] -= x;
  }
  template<class S>
  void operator+=(const S& x)
  {
    c[0] += x;
  }
  template<class S>
  void operator*=(const S& scale)
  {
    for (int i=0;i<POW2(Nvar);i++)
      c[i] *= scale;
  }
  template<class S>
  void operator/=(const S& scale)
  {
    for (int i=0;i<POW2(Nvar);i++)
      c[i] /= scale;
  }
  /*
  void operator/=(const ctaylor<T,Nvar>& t)
  {
    taylor<T,Nvar,Ndeg> tinv = 1/t;
    *this*=tinv;
  }
  */
  void operator*=(int scale)
  {
    for (int i=0;i<POW2(Nvar);i++)
      c[i] *= scale;
  }
  void operator*=(const ctaylor<T, Nvar>& t)
  {
    ctaylor_rec<T,Nvar>::multo(c,t.c);
  }
  /* Put sum_i coeff[i]*(this - this[0])^i in res,
     used when evaluating analytical functions of this */
  void compose(ctaylor<T, Nvar>& res, const T coeff[]) const
  {
    res = coeff[Nvar];
    for (int i=Nvar-1;i>=0;i--)
      {
	ctaylor_rec<T,Nvar>::multo_skipconst(res.c,c);
	res[0] += coeff[i];
      }
  }
};

template<class T, int Nvar>
static ctaylor<T, Nvar> operator*(const ctaylor<T, Nvar>& t1, 
				  const ctaylor<T, Nvar>& t2)
{
  ctaylor<T, Nvar> tmp = 0;
  ctaylor_rec<T,Nvar>::mul(tmp.c,t1.c,t2.c);
  return tmp;
}


// <> comparisons are taken to mean comparing the constant
// coefficient. This makes the transition from numbers to
// taylor objects easier.
template<class S, class T, int Nvar>
static bool operator<(const S &x, const ctaylor<T, Nvar> &t)
{
  return x < t[0];
}

template<class S, class T, int Nvar>
static bool operator>(const S &x, const ctaylor<T, Nvar> &t)
{
  return x > t[0];
}

template<class S, class T, int Nvar>
static bool operator<(const ctaylor<T, Nvar> &t, const S &x)
{
  return t[0] < x;
}

template<class S, class T, int Nvar>
static bool operator>(const ctaylor<T, Nvar> &t, const S &x)
{
  return t[0] > x;
}

template<class S, class T, int Nvar>
static bool operator!=(const ctaylor<T, Nvar> &t, const S &x)
{
  return t[0] != x;
}

template<class T, int Nvar, class S>
static ctaylor<T, Nvar> operator*(const S& x, const ctaylor<T, Nvar>& t)
{
  ctaylor<T, Nvar> tmp;
  for (int i=0;i<POW2(Nvar);i++)
    tmp[i] = x*t[i];
  return tmp;
}

template<class T, int Nvar, class S>
static ctaylor<T, Nvar> operator*(const ctaylor<T, Nvar>& t, const S& x)
{
  ctaylor<T, Nvar> tmp;
  for (int i=0;i<POW2(Nvar);i++)
    tmp[i] = x*t[i];
  return tmp;
}

template<class T, int Nvar, class S>
static ctaylor<T, Nvar> operator+(const S& x, const ctaylor<T, Nvar>& t)
{
  ctaylor<T, Nvar> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar, class S>
static ctaylor<T, Nvar> operator+(const ctaylor<T, Nvar>& t, const S& x)
{
  ctaylor<T, Nvar> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar>
static ctaylor<T, Nvar> operator+(const ctaylor<T, Nvar>& t1, 
				 const ctaylor<T, Nvar>& t2)
{
  ctaylor<T, Nvar> tmp;
  for (int i=0;i<POW2(Nvar);i++)
    tmp[i] = t1[i]+t2[i];
  return tmp;
}

template<class T, int Nvar>
static ctaylor<T, Nvar> operator-(const T& x, const ctaylor<T, Nvar>& t)
{
  ctaylor<T, Nvar> tmp = -t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar>
static ctaylor<T, Nvar> operator-(const ctaylor<T, Nvar>& t, const T &x)
{
  ctaylor<T, Nvar> tmp = t;
  tmp[0] -= x;
  return tmp;
}

template<class T, int Nvar>
static ctaylor<T, Nvar> operator-(int x, const ctaylor<T, Nvar>& t)
{
  ctaylor<T, Nvar> tmp = -t;
  tmp[0] += x;
  return tmp;
}

template<class T, int Nvar>
static ctaylor<T, Nvar> operator-(const ctaylor<T, Nvar>& t, int x)
{
  ctaylor<T, Nvar> tmp = t;
  tmp[0] -= x;
  return tmp;
}

template<class T, int Nvar>
static ctaylor<T, Nvar> operator-(const ctaylor<T, Nvar>& t1, 
				  const ctaylor<T, Nvar>& t2)
{
  ctaylor<T, Nvar> tmp;
  for (int i=0;i<POW2(Nvar);i++)
    tmp[i] = t1[i]-t2[i];
  return tmp;
}

#include "ctaylor_math.h"

#ifdef TAYLOR_CXXIO

#include <iostream>

template<class num, int Nvar>
static std::ostream &operator<<(std::ostream& stream, 
				const ctaylor<num,Nvar> &t)
{
  stream << "{" << t[0];
  for (int i=1;i<POW2(Nvar);i++)
    stream << ", " << t[i];
  stream << "}";
  return stream;
}

#endif

#endif
