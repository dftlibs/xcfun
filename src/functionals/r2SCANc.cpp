#include "functional.hpp"
#include "constants.hpp"
#include "SCAN_like_eps.hpp"

template<class num>
static num r2SCANC(const densvars<num> &d)
{
  num eps_c = SCAN_eps::SCAN_C(d, 2, 1, 1); 
 
  return eps_c;
}

FUNCTIONAL(XC_R2SCANC) = {
  "r2SCAN correlation functional",
  "r2SCAN correlation functional.\n"
  "Restored-Regularised SCAN functional\n"
  "J Furness, in preparation (arXiv:2008.03374)"
  "Implemented by James Furness\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(r2SCANC)
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
  XC_PARTIAL_DERIVATIVES,
};
