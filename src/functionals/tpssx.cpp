#include "functional.h"
#include "constants.h"
#include "tpssx_eps.h"

template<class num>
static num energy(const densvars<num> &d)
{
  num eps = tpssx_eps::tpssx_eps(d);
  return d.n*eps;
}

void setup_tpssx(functional &f)
{
  f.describe(XC_TPSSX, XC_MGGA,
	     "TPSS original exchange functional",
	     "TPSS original exchange functional.\n"
	     "J. Tao, J.P. Perdew, V. N. Staroverov, G. E. Scuseria,\n"
             "Climbing the Density Functional Ladder:\n"
	     "Nonempirical Meta-Generalized Gradient Approximation\n" 
             "Designed for Molecules and Solids,\n"
	     "Phys. Rev. Lett. 91 (2003) 146401\n"
	     "Implemented by Andrea Debnarova\n");
  SET_MGGA_ENERGY_FUNCTION(f,energy);
  // Test case from ??
 // const double d[1] = 
 //  {0};
 // const double ref[1] =
  //  {0};
 // f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}

