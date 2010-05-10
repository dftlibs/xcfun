#include "functional.h"
#include "pbex.h"

template<class num>
static num energy(const densvars<num> &d)
{
  return pbex::energy_pbe_ab(pbex::R_revpbe,d.a,d.gaa)                                              
    + pbex::energy_pbe_ab(pbex::R_revpbe,d.b,d.gbb);
}


void setup_revpbex(functional &f)
{
  f.describe("revpbex",XC_GGA,
	     "Revised PBE Exchange Functional",
	     "Revised PBE Exchange Functional\n"
	     "Y. Zhang and W., Phys. Rev. Lett 80, 890 (1998)\n"
	     "Implemented by Ulf Ekstrom and Andre Gomes\n");
  SET_GGA_ENERGY_FUNCTION(f,energy);
}
