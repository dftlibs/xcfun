#include "functional.hpp"
#include "constants.hpp"

//  von Weizsacker kinetic energy functional

template<class num>
static num tw(const densvars<num> &d)
{
    using xc_constants::CF;

    return 1./8.*pow(d.gaa+d.gbb,2.0)/d.n;
}

FUNCTIONAL(XC_TW) = {
  "von Weizsacker Kinetic Energy Functional",
  "von Weizsacker Kinetic Energy Functional\n"
  "Implemented by AB and SR.\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(tw)
  XC_A_B_GAA_GAB_GBB,
};

