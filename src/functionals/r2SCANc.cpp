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
  "J Furness, A Kaplan, J Ning, J Perdew, & J Sun; J. Chem. Phys. Lett.; Accepted (DOI: 10.1021/acs.jpclett.0c02405)"
  "Implemented by James Furness\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(r2SCANC)
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
  XC_PARTIAL_DERIVATIVES,
  1,
  1e-11,
  {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
  {-7.52536253477e-03, -1.99923310008e-02, -6.45082284064e-02,
   1.10376734989e-02, 2.20753469977e-02, 1.10376734989e-02,
   -1.76333423096e-02, -1.76333423096e-02}
};
