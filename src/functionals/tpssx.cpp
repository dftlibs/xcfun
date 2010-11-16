#include "functional.h"
#include "constants.h"
#include "tpssx_eps.h"

using namespace tpssx_eps;

template<class num>
static num ENERGY_FUNCTION(XC_TPSSX)(const densvars<num> &d)
{
  num Fxa = F_x(2*d.a,4*d.gaa,2*d.taua);
  num epsxunif_a = fx_unif(2*d.a);
  num Fxb = F_x(2*d.b,4*d.gbb,2*d.taub);
  num epsxunif_b = fx_unif(2*d.b);
  return (epsxunif_a*Fxa + epsxunif_b*Fxb)/2;
}

NEW_TMGGA_FUNCTIONAL(XC_TPSSX);
SHORT_DESCRIPTION(XC_TPSSX) = "TPSS original exchange functional";
LONG_DESCRIPTION(XC_TPSSX) =
	     "TPSS original exchange functional.\n"
	     "J. Tao, J.P. Perdew, V. N. Staroverov, G. E. Scuseria,\n"
             "Climbing the Density Functional Ladder:\n"
	     "Nonempirical Meta-Generalized Gradient Approximation\n" 
             "Designed for Molecules and Solids,\n"
	     "Phys. Rev. Lett. 91 (2003) 146401\n"
	     "Implemented by Andrea Debnarova\n";
TEST_VARS(XC_TPSSX) = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
TEST_MODE(XC_TPSSX) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_TPSSX) = 1;
TEST_THRESHOLD(XC_TPSSX) = 1e-8;
TEST_IN(XC_TPSSX) = 
    { 0.153652558932587,
      0.153652558932587,
      8.390981882024583E-002,
      8.390981882024583E-002,
      8.390981882024583E-002,
      6.826262722466284E-002,
      6.826262722466284E-002 };
TEST_OUT(XC_TPSSX) =
    {-1.7369862616505458e-01, //XCFun energy from ADF
     -0.696647334686932,
     -0.696647334686932,
     -0.104380091839662,
     0.0,
     -0.104380091839662,     
     0.128315600541666,       
     0.128315600541666};

