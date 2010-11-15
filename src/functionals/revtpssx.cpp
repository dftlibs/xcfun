#include "functional.h"
#include "constants.h"
#include "revtpssx_eps.h"

template<class num>
static num ENERGY_FUNCTION(XC_REVTPSSX)(const densvars<num> &d)
{
  return  revtpssx_eps::revtpssx_eps(d);
}      

NEW_TMGGA_FUNCTIONAL(XC_REVTPSSX);
SHORT_DESCRIPTION(XC_REVTPSSX) = "Reviewed TPSS exchange functional";
LONG_DESCRIPTION(XC_REVTPSSX) =
	     "Reviewed TPSS exchange functional.\n"
	     "J.P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, J. Sun,\n"
             "Workhorse Semilocal Density Functional\n"
	     "for Condensed Matter Physics and Quantum Chemistry\n" 
             "Phys. Rev. Lett. 103 (2009) 026403\n"
	     "Implemented by Andrea Debnarova\n";
NO_TEST(XC_REVTPSSX);

