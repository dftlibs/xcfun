#include "functional.h"
#include "pbex.h"

template<class num>
static num energy(const densvars<num> &d)
{
  using namespace pbex;
  using pw91_like_x_internal::prefactor;
  return prefactor(d.a)*enhancement_RPBE(d.a,d.gaa) +
    prefactor(d.b)*enhancement_RPBE(d.b,d.gbb);
}


void setup_rpbex(functional &f)
{
  f.describe(XC_RPBEX, XC_GGA,
	     "RPBE Exchange Functional",
	     "RPBE Exchange Functional\n"
	     "Hammer, B. Hansen, L.B., Norskov, J.K.; PRB (59) p.7413, 1999\n"
	     "Implemented by Ulf Ekstrom\n");
  SET_GGA_ENERGY_FUNCTION(f,energy);
}
