#include "functional.h"
#include <iostream>
using namespace std;

template<class num>
static num becke_alpha(const num &na, const num &gaa)
{
  static const parameter c = pow(81/(4*M_PI),1.0/3.0)/2;
  static const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num lda = -c*na43;
  num chi = sqrt(gaa)/na43;
  num b88 = -(d*na43*chi*chi)/(1+6*d*chi*asinh(chi));
  return lda + b88;
}

template<class num>
static num becke_corr(const num &na, const num &gaa)
{
  static const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num chi = sqrt(gaa)/na43;
  return -(d*na43*chi*chi)/(1+6*d*chi*asinh(chi));
}

// Short range becke exchange as used in camb3lyp If mu=0 this reduces
// to the standard beckex for which must be used instead.
// FIXME: The erf + something is basically the erf minus
// its asymptotic expansion. This is horribel for numerics,
// will have to code a special function. 
// As coded here the code will fail if mu = 0, in which case the
// regular beckex should be used.

template<class num>
static num becke_sr(parameter mu, const num &na, const num &gaa)
{
  static const parameter cparam = pow(81/(4*M_PI),1.0/3.0)/2;
  static const parameter d = 0.0042;
  num na43 = pow(na,4.0/3.0);
  num chi = sqrt(gaa)/na43;
  num K = 2*(cparam + (d*chi*chi)/(1+6*d*chi*asinh(chi)));
  num a = mu*sqrt(K)/(6*M_PI*pow(na,1.0/3.0));
  num b = exp(-1/(4*a*a))-1;
  num c = 2*a*a*b + 0.5;
  return -0.5*na43*K*(1-8.0/3.0*a*(sqrt(M_PI)*erf(1/(2*a))+2*a*(b-c)));
}

template<class num>
static num energy(const densvars<num> &d)
{
  return becke_alpha(d.a,d.gaa) + becke_alpha(d.b,d.gbb);
}

template<class num>
static num energy_corr(const densvars<num> &d)
{
  return becke_corr(d.a,d.gaa) + becke_corr(d.b,d.gbb);
}

template<class num>
static num energy_sr(const densvars<num> &d)
{
  parameter mu = 1e-10;
  return becke_sr(mu,d.a,d.gaa) + becke_sr(mu,d.b,d.gbb);
}

void setup_beckex(functional &f)
{
  f.describe("beckex",XC_GGA,
	     "Becke 88 exchange",
	     "Becke 88 exchange including Slater part\n"
	     "A.D. Becke, Density-functional exchange-energy approximation\n"
	     "with correct asymptotic behaviour, Phys. Rev. A38 (1988) 3098-3100.\n"
	     "Implemented by Ulf Ekstrom\n"
	     "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n");

  SET_GGA_ENERGY_FUNCTION(f,energy);

  // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_b88.html
  static const double d[5] = 
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
  static const double ref[21] =
    {	-0.277987329958E+03,
	-0.385951846654E+01,
	-0.381309494319E+01,
	-0.172434478018E-04,
	-0.173712338362E-04,
	0.000000000000E+00,
	-0.441426807406E-01,
	0.000000000000E+00,
	0.201415922856E-06,
	0.000000000000E+00,
	0.000000000000E+00,
	-0.447245742260E-01,
	0.000000000000E+00,
	0.195961359539E-06,
	0.000000000000E+00,
	0.700742719647E-11,
	0.000000000000E+00,
	0.000000000000E+00,
	0.718678968862E-11,
	0.000000000000E+00,
	0.000000000000E+00
    };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}

void setup_beckexcorr(functional &f)
{
  f.describe("beckexcorr",XC_GGA,
	     "Becke 88 exchange correction",
	     "Becke 88 exchange not including Slater part\n"
	     "A.D. Becke, Density-functional exchange-energy approximation\n"
	     "with correct asymptotic behaviour, Phys. Rev. A38 (1988) 3098-3100.\n"
	     "Implemented by Ulf Ekstrom\n"
	     "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n");

  SET_GGA_ENERGY_FUNCTION(f,energy_corr);
}

void setup_beckexsr(functional &f)
{
  f.describe("beckex_sr",XC_GGA,
	     "Short range Becke 88 exchange",
	     "Implemented by Ulf Ekstrom\n"
	     "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n");

  SET_GGA_ENERGY_FUNCTION(f,energy_sr);

}
