#include "functional.hpp"
#include "constants.hpp"
#include "SCAN_like_eps.hpp"

template<class num>
static num rppSCANC(const densvars<num> &d)
{
  num eps_c = SCAN_eps::SCAN_C(d, 2, 1, 0); 
 
  return eps_c;
}

FUNCTIONAL(XC_RPPSCANC) = {
  "r++SCAN correlation functional",
  "r++SCAN correlation functional.\n"
  "Restored-Regularised SCAN functional\n"
  "J Furness, in preparation"
  "Implemented by James Furness\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(rppSCANC)
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
  XC_PARTIAL_DERIVATIVES,
};
