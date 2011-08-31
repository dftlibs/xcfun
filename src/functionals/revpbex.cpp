#include "functional.hpp"
#include "pbex.hpp"

template<class num>
static num revpbex(const densvars<num> &d)
{
  return pbex::energy_pbe_ab(pbex::R_revpbe,d.a,d.gaa)
    + pbex::energy_pbe_ab(pbex::R_revpbe,d.b,d.gbb);
}

FUNCTIONAL(XC_REVPBEX) = {
  "Revised PBE Exchange Functional",
  "Revised PBE Exchange Functional\n"
  "Y. Zhang and W., Phys. Rev. Lett 80, 890 (1998)\n"
  "Implemented by Ulf Ekstrom and Andre Gomes\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(revpbex)
};

