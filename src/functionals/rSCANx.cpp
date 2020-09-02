#include "functional.hpp"
#include "constants.hpp"
#include "SCAN_like_eps.hpp"

template<class num>
static num rSCANX(const densvars<num> &d)
{
  num Fx_a;
  num Fx_b;
  Fx_a = SCAN_eps::get_SCAN_Fx(2.0*d.a, 4.0*d.gaa, 2.0*d.taua, 1, 1, 0); 
  num epsxunif_a = SCAN_eps::fx_unif(2*d.a);
  Fx_b = SCAN_eps::get_SCAN_Fx(2.0*d.b, 4.0*d.gbb, 2.0*d.taub, 1, 1, 0); 
  num epsxunif_b = SCAN_eps::fx_unif(2*d.b);

  return 0.5*(Fx_a*epsxunif_a + Fx_b*epsxunif_b);
}

FUNCTIONAL(XC_RSCANX) = {
  "rSCAN exchange functional",
  "rSCAN exchange functional.\n"
  "Regularised SCAN functional\n"
  "[1] A. P. Bartók and J. R. Yates, J. Chem. Phys. 150, 161101 (2019)."
  "Implemented by James Furness\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(rSCANX)
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
  XC_PARTIAL_DERIVATIVES,
};
