#include "functional.hpp"
#include "pbex.hpp"

template<class num>
static num rpbex(const densvars<num> &d)
{
  using namespace pbex;
  using pw91_like_x_internal::prefactor;
  return prefactor(d.a)*enhancement_RPBE(d.a,d.gaa) +
    prefactor(d.b)*enhancement_RPBE(d.b,d.gbb);
}


FUNCTIONAL(XC_RPBEX) = {
  "RPBE Exchange Functional",
  "RPBE Exchange Functional\n"
  "Hammer, B. Hansen, L.B., Norskov, J.K.; PRB (59) p.7413, 1999\n"
  "Implemented by Ulf Ekstrom\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(rpbex)
};


