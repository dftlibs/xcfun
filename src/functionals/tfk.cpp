#include "functional.h"
#include "constants.h"

//  Thomas-Fermi kinetic energy functional

template<class num>
static num tfk(const densvars<num> &d)
{
    using xc_constants::CF;

    return CF*pow(d.n, 5.0/3.0);
}

FUNCTIONAL(XC_TFK) = {
  "Thomas-Fermi Kinetic Energy Functional",
  "Thomas-Fermi Kinetic Energy Functional\n"
  "Implemented by Andre Gomes.\n",
  XC_DENSITY,
  ENERGY_FUNCTION(tfk)
  XC_A_B,
  XC_PARTIAL_DERIVATIVES,
  1,
  1e-5,
  {1., .8},
  { 7.64755771625168, 7.08107195949229, 7.08107195949229, 2.62261924425641, 2.62261924425641, 2.62261924425641 }
};

