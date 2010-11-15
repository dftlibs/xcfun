#include "functional.h"
#include "pbex.h"

template<class num>
static num ENERGY_FUNCTION(XC_RPBEX)(const densvars<num> &d)
{
  using namespace pbex;
  using pw91_like_x_internal::prefactor;
  return prefactor(d.a)*enhancement_RPBE(d.a,d.gaa) +
    prefactor(d.b)*enhancement_RPBE(d.b,d.gbb);
}


NEW_GGA_FUNCTIONAL(XC_RPBEX);
SHORT_DESCRIPTION(XC_RPBEX) = "RPBE Exchange Functional";
LONG_DESCRIPTION(XC_RPBEX) =	     "RPBE Exchange Functional\n"
	     "Hammer, B. Hansen, L.B., Norskov, J.K.; PRB (59) p.7413, 1999\n"
	     "Implemented by Ulf Ekstrom\n";
NO_TEST(XC_RPBEX);

