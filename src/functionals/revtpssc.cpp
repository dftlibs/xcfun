#include "functional.hpp"
#include "constants.hpp"
#include "revtpssc_eps.hpp"

template<class num>
static num revtpssc(const densvars<num> &d)
{
  num eps = revtpssc_eps::revtpssc_eps(d);
  return d.n*eps;
}

FUNCTIONAL(XC_REVTPSSC) = {
  "Revised TPSS correlation functional",
  "Revised TPSS correlation functional.\n"
  "J.P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, J. Sun,\n"
  "Workhorse Semilocal Density Functional\n"
  "for Condensed Matter Physics and Quantum Chemistry\n" 
  "Phys. Rev. Lett. 103 (2009) 026403\n"
  "Implemented by Andrea Debnarova\n",
  XC_DENSITY | XC_GRADIENT | XC_KINETIC,
  ENERGY_FUNCTION(revtpssc)
};


