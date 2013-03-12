#include "functional.hpp"

template<class num>
static num vW_alpha(const num &na, const num &gaa)
{
  return gaa/(8*na);
}

template<class num>
static num vW(const densvars<num> &d)
{
  return vW_alpha(d.a,d.gaa) + vW_alpha(d.b,d.gbb);
}

FUNCTIONAL(XC_VWK) = {
  "von Weizsaecker kinetic energy",
  "von Weizsaecker kinetic energy\n"
  "Implemented by Borgoo/Ekstrom.\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(vW)
};
 
