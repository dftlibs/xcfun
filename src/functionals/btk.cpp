#include "functional.hpp"

template<class num>
static num btk_alpha(const num &na, const num &gaa)
{
  const parameter qav = 0.3434125;
  const parameter beta = 1.990328;
  const parameter fudge = 1e-24;
  num na53 = pow(na,5.0/3.0);
  return beta*na53*pow(fudge+gaa,0.5*qav)/pow(na,4.0/3.0*qav);
}

template<class num>
static num btk(const densvars<num> &d)
{
  return 0.5*(btk_alpha(2*d.a,4*d.gaa) + btk_alpha(2*d.b,4*d.gbb));
}

FUNCTIONAL(XC_BTK) = {
  "Borgoo-Tozer TS",
  "Borgoo-Tozer kinetic energy functional\n"
  "Implemented by Borgoo/Ekstrom.\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(btk)
};
 
