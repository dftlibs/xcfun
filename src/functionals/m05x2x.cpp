#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M05-2X exchange functional. to be used with HF exchange factor of .56 

template<class num>
static num ENERGY_FUNCTION(XC_M05X2X) (const densvars<num> &d)
{
   using m0xy_metagga_xc_internal::fw;

   const parameter param_a[12] =
     {  1.000000e+00, -5.683300e-01, -1.300570e+00,  5.500700e+00,  9.064020e+00,
       -3.221075e+01, -2.373298e+01,  7.022996e+01,  2.988614e+01, -6.025778e+01,
       -1.322205e+01,  1.523694e+01 };

   return (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua)
           + pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) );
}

NEW_TMGGA_FUNCTIONAL(XC_M05X2X);
SHORT_DESCRIPTION(XC_M05X2X) = "M05-2X exchange";
LONG_DESCRIPTION(XC_M05X2X) =             "M05-2X Meta-Hybrid Exchange Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n";
TEST_VARS(XC_M05X2X) = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
TEST_MODE(XC_M05X2X) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_M05X2X) = 1;
TEST_THRESHOLD(XC_M05X2X) = 3e-5;
TEST_IN(XC_M05X2X) = {1., .8, 1., 1., 1.,  0.165,   0.1050};
TEST_OUT(XC_M05X2X) =    { -1.38233309, -0.19638222, -0.08614105, -0.00289174,   
      0.00000000, -0.00365982, -3.18842316000000,  -3.90587738 }; 

