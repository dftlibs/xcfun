#include "functional.h"
#include "constants.h"
#include "tpssc_eps.h"

template<class num>
static num ENERGY_FUNCTION(XC_TPSSC)(const densvars<num> &d)
{
  num eps = tpssc_eps::tpssc_eps(d);
  return d.n*eps;
} 

NEW_TMGGA_FUNCTIONAL(XC_TPSSC);
SHORT_DESCRIPTION(XC_TPSSC) ="TPSS original correlation functional";
LONG_DESCRIPTION(XC_TPSSC) =
	     "TPSS original correlation functional.\n"
	     "J. Tao, J.P. Perdew, V. N. Staroverov, G. E. Scuseria,\n"
             "Climbing the Density Functional Ladder:\n"
	     "Nonempirical Meta-Generalized Gradient Approximation\n" 
             "Designed for Molecules and Solids,\n"
	     "Phys. Rev. Lett. 91 (2003) 146401\n"
	     "Implemented by Andrea Debnarova\n";
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
TEST_VARS(XC_TPSSC) = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
TEST_MODE(XC_TPSSC) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_TPSSC) = 1;
TEST_THRESHOLD(XC_TPSSC) = 1e-6;
TEST_IN(XC_TPSSC) = {1,2,3,4,5,6,7};
TEST_OUT(XC_TPSSC) = {-2.1824017471364521e-01,
			-0.114815778042036,
			-7.968561473875205E-002,
			8.304923723228601E-004,
			1.671559231642408E-003,
			8.327578102972385E-004,
			-1.440645713140351E-005,
			-1.440645713140351E-005};

