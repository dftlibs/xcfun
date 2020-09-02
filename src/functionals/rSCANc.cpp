#include "functional.hpp"
#include "constants.hpp"
#include "SCAN_like_eps.hpp"

template<class num>
static num rSCANC(const densvars<num> &d)
{
  num eps_c = SCAN_eps::SCAN_C(d, 1, 1, 0); 
 
  return eps_c;
}

FUNCTIONAL(XC_RSCANC) = {
  "rSCAN correlation functional",
  "rSCAN correlation functional.\n"
  "Regularised SCAN functional\n"
  "A. P. Bartók and J. R. Yates, J. Chem. Phys. 150, 161101 (2019)."
  "Implemented by James Furness\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(rSCANC)
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
  XC_PARTIAL_DERIVATIVES,
};
