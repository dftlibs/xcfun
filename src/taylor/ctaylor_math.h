#include "taylor.h"

template<class T,int Nvar, class S>
static ctaylor<T,Nvar> operator/(const S &x, const ctaylor<T,Nvar>& t)
{
  taylor<T,1,Nvar> tmp;
  inv_taylor(tmp,t[0]);
  tmp*=x;
  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}

template<class T,int Nvar>
static ctaylor<T,Nvar> operator/(const ctaylor<T,Nvar>&t1, 
				 const ctaylor<T,Nvar>&t2)
{
  taylor<T,1,Nvar> tmp;
  inv_taylor(tmp,t2[0]);
  ctaylor<T,Nvar> res;
  t2.compose(res,tmp.c);
  res*=t1;
  return res;
}

template<class T,int Nvar>
static ctaylor<T,Nvar> operator/(const ctaylor<T,Nvar>& t, 
				const T &x)
{
  ctaylor<T,Nvar> tmp = t;
  tmp *= 1/x;
  return tmp;
}

template<class T, int Nvar, class S>
static ctaylor<T,Nvar> operator/(const ctaylor<T,Nvar>& t, const S &x)
{
  ctaylor<T,Nvar> tmp = t;
  tmp *= 1/T(x);
  return tmp;
}

template<class T,int Nvar>
static ctaylor<T,Nvar> exp(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  exp_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}


// exp(x)-1, but accurate for small x
template<class T,int Nvar>
static ctaylor<T,Nvar> expm1(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  exp_taylor(tmp,t[0]);

  // Only constant value is affected by the cancellation
  if (fabs(t[0]) > 1e-3)
    tmp[0] -= 1;
  else
    tmp[0] = 2*exp(t[0]/2)*sinh(t[0]/2);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}



template<class T,int Nvar>
static ctaylor<T,Nvar> log(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  log_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}

// We need this version with double a argument to prevent truncation
// to int.
template<class T,int Nvar>
static ctaylor<T,Nvar> pow(const ctaylor<T,Nvar> &t, const double &a)
{
  taylor<T,1,Nvar> tmp;
  pow_taylor(tmp,t[0],a);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}


template<class T,int Nvar>
static ctaylor<T,Nvar> sqrt(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  sqrt_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}


template<class T,int Nvar>
static ctaylor<T,Nvar> cbrt(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  cbrt_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}


// Integer exponent version is analytical at t[0] = 0
// This function gets priority over the normal pow 
// when the exponent is an integer, but does not force
// conversion to integer.
template<class T,int Nvar>
static ctaylor<T,Nvar> pow(const ctaylor<T,Nvar> &t, int n)
{
  if (n > 0)
    {
      ctaylor<T,Nvar> res = t;
      while (n-- > 1)
	res *= t;
      return res;
    }
  else if (n < 0)
    {
      return pow(t,double(n));
    }
  else 
    {
      ctaylor<T,Nvar> res(1);
      return res;
    }
}

template<class T,int Nvar>
static ctaylor<T,Nvar> atan(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  atan_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}


template<class T,int Nvar>
static ctaylor<T,Nvar> erf(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  erf_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}



template<class T,int Nvar>
static ctaylor<T,Nvar> sin(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  sin_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}



template<class T,int Nvar>
static ctaylor<T,Nvar> cos(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  cos_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}



template<class T,int Nvar>
static ctaylor<T,Nvar> asinh(const ctaylor<T,Nvar> &t)
{
  taylor<T,1,Nvar> tmp;
  asinh_taylor(tmp,t[0]);

  ctaylor<T,Nvar> res;
  t.compose(res,tmp.c);
  return res;
}



/*
  The original function is unstable for small t[0] values similarly to
  the Boys function. Use an [8,8] Pade approximation when |t[0]| is
  small. This works less well but still ok in single precision.
 */
template<int Nvar>
static ctaylor<double,Nvar> sqrtx_asinh_sqrtx(const ctaylor<double,Nvar> &t)
{
  assert(t[0] > -0.5);
  // Coefficients of an [8,8] Pade approximation at x = 0
  static const double P[] = {0,
			     3.510921856028398e3,
			     1.23624388373212e4,
			     1.734847003883674e4,
			     1.235072285222234e4,
			     4.691117148130619e3,
			     9.119186273274577e2,
			     7.815848629220836e1,	     
			     1.96088643023654e0};
  static const double Q[] = {3.510921856028398e3,
			     1.29475924799926e4,
			     1.924308297963337e4,
			     1.474357149568687e4,    
			     6.176496729255528e3,
			     1.379806958043824e3,
			     1.471833349002349e2,
			     5.666278232986776e0,
			     2.865104054302032e-2};
  if (fabs(t[0]) < 0.5)
    {
      // Shift polys
      assert(Nvar<=8);
      taylor<double,1,Nvar> pp,pq,ppade;
      reinterpret_cast<const taylor<double,1,8> *>(P)->shift(pp,&t[0]);
      reinterpret_cast<const taylor<double,1,8> *>(Q)->shift(pq,&t[0]);
      ppade = pp/pq;
      ctaylor<double,Nvar> res;
      t.compose(res,ppade.c);
      return res;
    }
    else
    {
      // This is the unstable form
      ctaylor<double,Nvar> s = sqrt(t);
      return s*asinh(s);
    }
}


template<class T,int Nvar>
static ctaylor<T,Nvar> min(const ctaylor<T,Nvar> &a,
			   const ctaylor<T,Nvar> &b)
{
  if (a <= b)
    return a;
  else
    return b;
}

template<class T,int Nvar>
static ctaylor<T,Nvar> max(const ctaylor<T,Nvar> &a,
			   const ctaylor<T,Nvar> &b)
{
  if (a > b)
    return a;
  else
    return b;
}


