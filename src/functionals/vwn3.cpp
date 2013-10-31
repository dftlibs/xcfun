#include "functional.hpp"
#include "vwn.hpp"

template<class num>
static num vwn3c(const densvars<num> &d)
{
  return d.n*vwn::vwn3_eps(d);
}

FUNCTIONAL(XC_VWN3C) = 
{
  "VWN3 LDA Correlation functional",
  "VWN3 LDA Correlation functional\n"
  "S.H. Vosko, L. Wilk, and M. Nusair: Accurate spin-dependent\n"
  "electron liquid correlation energies for local spin density\n"
  "calculations: a critical analysis, Can. J. Phys. 58 (1980) 1200-1211.\n"
  "Originally from Dalton, polished and converted by Ulf Ekstrom.\n",
  XC_DENSITY,
  ENERGY_FUNCTION(vwn3c)
};

