#include "functional.h"
#include "constants.h"
#include "revtpssc_eps.h"

template<class num>
static num ENERGY_FUNCTION(XC_REVTPSSC)(const densvars<num> &d)
{
  num eps = revtpssc_eps::revtpssc_eps(d);
  return d.n*eps;
}

NEW_TMGGA_FUNCTIONAL(XC_REVTPSSC);
SHORT_DESCRIPTION(XC_REVTPSSC) = "Revised TPSS correlation functional";
LONG_DESCRIPTION(XC_REVTPSSC) =
	     "Revised TPSS correlation functional.\n"
	     "J.P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, J. Sun,\n"
             "Workhorse Semilocal Density Functional\n"
	     "for Condensed Matter Physics and Quantum Chemistry\n" 
             "Phys. Rev. Lett. 103 (2009) 026403\n"
             "Implemented by Andrea Debnarova\n";
NO_TEST(XC_REVTPSSC);

