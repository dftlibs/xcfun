#include "functional.h"
#include "constants.h"
#include "revtpssx_eps.h"

template<class num>
static num energy(const densvars<num> &d)
{
  num eps = revtpssx_eps::revtpssx_eps(d);
  return eps; 
}      

void setup_revtpssx(functional &f)
{
  f.describe(XC_REVTPSSX, XC_MGGA,
	     "Reviewed TPSS exchange functional",
	     "Reviewed TPSS exchange functional.\n"
	     "J.P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, J. Sun,\n"
             "Workhorse Semilocal Density Functional\n"
	     "for Condensed Matter Physics and Quantum Chemistry\n" 
             "Phys. Rev. Lett. 103 (2009) 026403\n"
	     "Implemented by Andrea Debnarova\n");
  SET_MGGA_ENERGY_FUNCTION(f,energy);
  // Test case from ??
 // const double d[1] = 
 //  {0};
 // const double ref[1] =
  //  {0};
 // f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}

