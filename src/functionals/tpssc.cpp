#include "functional.h"
#include "constants.h"
#include "tpssc_eps.h"

template<class num>
static num energy(const densvars<num> &d)
{
  num eps = tpssc_eps::tpssc_eps(d);
  return d.n*eps;
} 

void setup_tpssc(functional &f)
{
  f.describe(XC_TPSSC, XC_MGGA,
	     "TPSS original correlation functional",
	     "TPSS original correlation functional.\n"
	     "J. Tao, J.P. Perdew, V. N. Staroverov, G. E. Scuseria,\n"
             "Climbing the Density Functional Ladder:\n"
	     "Nonempirical Meta-Generalized Gradient Approximation\n" 
             "Designed for Molecules and Solids,\n"
	     "Phys. Rev. Lett. 91 (2003) 146401\n"
	     "Implemented by Andrea Debnarova\n");
  SET_MGGA_ENERGY_FUNCTION(f,energy);
  /* // Test case from A**
  const double d[] = {0.153652558932587,
		      0.153652558932587,
		      8.390981882024583E-002,
		      8.390981882024583E-002,
		      8.390981882024583E-002,
		      6.826262722466284E-002,
		      6.826262722466284E-002};
  const double ref[] = {0,//Fill in energy
			   -7.299596527683171E-002,
		       -7.299596527683171E-002,
		       1.247652223035057E-002, // Probably incorrect
		       5.978495466748646E-002, // Probably incorrect
		       1.247652223035057E-002, // Probably incorrect
		       -4.068185412322921E-002,
		       -4.068185412322921E-002};
*/
  const double d2[] = {1,2,3,4,5,6,7};
  const double ref2[] = {-2.1824017471364521e-01,
			-0.114815778042036,
			-7.968561473875205E-002,
			8.304923723228601E-004,
			1.671559231642408E-003,
			8.327578102972385E-004,
			-1.440645713140351E-005,
			-1.440645713140351E-005};
 f.add_test(XC_VARS_AB,1,d2,ref2,1e-6);
}

