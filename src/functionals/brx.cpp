#ifdef TODO
/* Needs some work to work with ctaylor */
#include "functional.hpp"
#include "slater.hpp"

// This is the function we want to find roots for
template<class T>
T BR_y(const T &x)
{
  return x*exp(-2.0/3.0*x)/(x-2);
}


// Determine x as a function of y = BR_y(x), using first 
// a guess and then root polishing.
// We need about 3-5 NR iterations depending on the y value.
ireal_t BR(ireal_t y)
{
  ireal_t d;
  taylor<ireal_t,1,2> x(0,0),fx;
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
      // d = (y-fx[0])/fx[1];
      // Halley:
      ireal_t yn = fx[0] - y;
      d = -yn*fx[1]/(fx[1]*fx[1] - yn*fx[2]);
      x[0] += d;
      if (++niter > 100)
	{
	  assert(0 && "Newton not converging in BR"); 
	  return 0;
	}
    }
  while (fabs(d)>1e-12);
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

template<class num>
static num polarized(const num &na, 
		     const num &gaa,
		     const num &lapa, 
		     const num &taua) // Becke tau here, no factor 1/2
{
  const parameter gam = 0.80;
  num D = taua - 0.25*gaa/na;
  num Q = 1.0/6.0*(lapa + 2*gam*D);
  num x = BR(2.0/3.0*pow(M_PI,2.0/3.0)*pow(na,5.0/3.0)/Q);
  num b = cbrt(pow3(x)*exp(-x)/(8*M_PI*na));
  return -(1-(1+0.5*x)*exp(-x))/b; //FIXME: use expm1
}

template<class num>
static num ENERGY_FUNCTION(XC_BRX)(const densvars<num> &d)
{
  return 0.5*(d.a*polarized(d.a,d.gaa,d.lapa,2*d.taua) +
	      d.a*polarized(d.b,d.gbb,d.lapb,2*d.taub));
}


NEW_LTMGGA_FUNCTIONAL(XC_BRX);
SHORT_DESCRIPTION(XC_BRX) = "BR exchange\n";
LONG_DESCRIPTION(XC_BRX) =
	     "Becke-Roussels exchange functional\n"
	     "Implemented by Ulf Ekstrom\n";
NO_TEST(XC_BRX);

#endif
