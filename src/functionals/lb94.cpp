#include "functional.h"
#include "vwn.h"
#include "slater.h"

//Have to be a bit special here since lb94 is not a functional but 
//a potential.

template<class T, int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> lb94(const densvars<taylor<T,Nvar,Ndeg> > &d)
{
  static const parameter beta = 0.05;
  taylor<T,Nvar,Ndeg> res = vwn::vwn5c(d) + slaterx(d);
  // LB94 is basically LDA with a GGA modification _for the potential_
  if (Ndeg >= 1)
    {
      // Alpha potential
      T x = d.gaa[0]*pow(d.a[0],-4.0/3.0);
      T v = -beta*pow(d.a[0],1.0/3.0)*x*x/(1+3*beta*x*asinh(x)); 
      // Use the chain rule to calculate contributions to
      // each dE/dx for the unknown x. 
      for (int i=0;i<Nvar;i++)
	res[1+i] += v*d.a[1+i];

      // Beta potential
      x = d.gbb[0]*pow(d.b[0],-4.0/3.0);
      v = -beta*pow(d.b[0],1.0/3.0)*x*x/(1+3*beta*x*asinh(x)); 
      for (int i=0;i<Nvar;i++)
	res[1+i] += v*d.b[1+i];
    }
  return res;
}

void setup_lb94(functional &f)
{
  f.describe("lb94",XC_GGA,
	     "LB94 Exchange-correlation functional",
	     "Exchange-correlation potential with correct asymtotic behavior\n"
	     "R. van Leeuwen and E. J. Baerends PRA 49, 2421 (1994)\n"
	     "Note that the LB94 energy is not well defined, here its"
	     "simply the SVWN5 energy, also used for order >= 2.\n"
	     "Implemented by Ulf Ekstrom.\n");
  SET_GGA_ENERGY_FUNCTION(f,lb94);
}
