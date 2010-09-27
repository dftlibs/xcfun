#include <iostream>
#include <cassert>
#include <cmath>
#define TAYLOR_CXXIO
#include "taylor.h"

/*
  Initial tests for the Becke Roussel functionals.
  This shows how to taylor expand the implicitly
  defined quantity x, see BR() below.

  by ulfek
 */

using namespace std;

// This is the function we want to find roots for
template<class T>
T BR_y(const T &x)
{
  return x*exp(-2.0/3.0*x)/(x-2);
}


// Determine x as a function of y = BR_y(x), using first 
// a guess and then root polishing.
// We need about 3-5 NR iterations depending on the y value.
double BR(double y)
{
  double d;
  taylor<double,1,2> x(0,0),fx;
  // More or less clever starting guesses
  if (y < 0)
    {
      if (y > -0.65)
	x[0] = -2*y;
      else
	x[0] = (6*exp(4.0/3.0)*y + 8)/(3*exp(4.0/3.0)*y+1);
    }
  else
    {
      if (y > 0.1)
	{
	  x[0] = (6*exp(4.0/3.0)*y + 8)/(3*exp(4.0/3.0)*y+1); 
	  if (y < 0.7)
	    x[0] -= 0.7*exp(-5*y); // Pragmatic correction
	}
      else
	{
	  x[0] = -3.0/2.0*log(y) + y*(47 + y*(-850 + y*4800));
	}
    }
  //Polish root with Newton-Raphson (or Halley's method)
  int niter = 0;
  do
    {
      fx = BR_y(x);
      // Newton:
       d = (y-fx[0])/fx[1];
      // Halley:
      double yn = fx[0] - y;
      //d = -yn*fx[1]/(fx[1]*fx[1] - yn*fx[2]);
      x[0] += d;
      if (++niter > 100)
	{
	  cerr << niter << " iterations reached, giving up. y = " << y << endl;
	  return 0;
	}
    }
  while (fabs(d)>1e-12);
  cerr << "niter = " << niter << endl;
  return x[0];
}


// Obtain the Taylor expansion of x(y), which is the
// inverse of BR_y. Use linear method for simplicity.
template<class T, int Ndeg>
void BR_taylor(const T &y0, taylor<T,1,Ndeg> &t)
{
  taylor<T,1,Ndeg> f,d;
  t = 0;
  t[0] = BR(y0);
  t[1] = 1;
  f = BR_y(t);
  t[1] = 1/f[1];
  // Linear method, for quadratic see i.e. Brent & Kung ~197x
  for (int i=2;i<=Ndeg;i++) 
    {
      f = BR_y(t);
      t[i] = -f[i]*t[1];
    }
}

/* This is a fully differentiable solver for Eq.(21) in 
   Becke and Roussel, PRA 39, 1989. t is the right hand
   side value, x is returned.
 */ 
template<class T,int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> BR(const taylor<T,Nvar,Ndeg> &t)
{
  taylor<T,1,Ndeg> tmp;
  BR_taylor(t[0],tmp);

  taylor<T,Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


int main()
{

  taylor<double,1,5> seed(3,0);

  for (double y = -100;y<100;y+= 5.001)
    cout << y << " " << BR(y) << endl;
  cout << "# taylor expansion coeffs at y = 3: " << BR(seed) << endl;
  return 0;
}
