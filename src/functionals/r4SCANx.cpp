#include "functional.hpp"
#include "constants.hpp"
#include "SCAN_like_eps.hpp"

template<class num>
static num r4SCANX(const densvars<num> &d)
{
  num Fx_a;
  num Fx_b;
  Fx_a = SCAN_eps::get_SCAN_Fx(2.0*d.a, 4.0*d.gaa, 2.0*d.taua, 2, 1, 2); 
  num epsxunif_a = SCAN_eps::fx_unif(2*d.a);
  Fx_b = SCAN_eps::get_SCAN_Fx(2.0*d.b, 4.0*d.gbb, 2.0*d.taub, 2, 1, 2); 
  num epsxunif_b = SCAN_eps::fx_unif(2*d.b);

  return 0.5*(Fx_a*epsxunif_a + Fx_b*epsxunif_b);
}

FUNCTIONAL(XC_R4SCANX) = {
  "r4SCAN exchange functional",
  "r4SCAN exchange functional.\n"
  "Restored-Regularised SCAN functional\n"
  "J Furness, in preparation"
  "Implemented by James Furness\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(r4SCANX)
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
  XC_PARTIAL_DERIVATIVES,
  1,
  1e-11,
  {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
  {-1.62076202198e-01, -8.60008367983e-01, -5.57023311096e-01, 
      -3.79265181913e-02, 0.00000000000e+00, -1.04606093722e-01,
      5.89149583569e-02, 5.50319523444e-02}
};
