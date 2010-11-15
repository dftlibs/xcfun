#include "functional.h"
#include "pbex.h"

template<class num>
static num ENERGY_FUNCTION(XC_REVPBEX)(const densvars<num> &d)
{
  return pbex::energy_pbe_ab(pbex::R_revpbe,d.a,d.gaa)
    + pbex::energy_pbe_ab(pbex::R_revpbe,d.b,d.gbb);
}


NEW_GGA_FUNCTIONAL(XC_REVPBEX);
SHORT_DESCRIPTION(XC_REVPBEX) = "Revised PBE Exchange Functional",;
LONG_DESCRIPTION(XC_REVPBEX) =
	     "Revised PBE Exchange Functional\n"
	     "Y. Zhang and W., Phys. Rev. Lett 80, 890 (1998)\n"
	     "Implemented by Ulf Ekstrom and Andre Gomes\n");
NO_TEST(XC_REVPBEX);

