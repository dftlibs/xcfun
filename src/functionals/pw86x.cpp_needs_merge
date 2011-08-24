#include "functional.h"

template<class num> 
static num pw86x(const num &na, const num &gaa)
{
  const parameter a = 1.0;
  const parameter b = 1.296;
  const parameter c = 14.0;
  const parameter d = 0.20;
  num rho = 2*na, grad2 = 4*gaa;
  const num Ax = -pow(3.0/M_PI,1.0/3.0)*3.0/4.0;
  const num kf = pow(3.0*pow(M_PI,2)*rho,1.0/3.0);
  num s2 = grad2/pow(2.0*kf*rho,2);
  num F = pow(a + s2*(b + s2*(c + d*s2)), 1.0/15.0);
  return Ax*pow(rho,4.0/3.0)*F;
}

template<class num>
static num energy(const densvars<num> &d)
{
  return 0.5*(pw86x(d.a,d.gaa) + pw86x(d.b,d.gbb));
}

void setup_pw86x(functional &f)
{
  f.describe(XC_PW86X, XC_GGA,
	     "Perdew-Wang 1986 GGA Exchange Functional",
	     "Perdew-Wang 1986 GGA Exchange Functional\n"
	     "Phys. Rev. B 33. 8800 (1986)"); 
  SET_GGA_ENERGY_FUNCTION(f,energy);
}
