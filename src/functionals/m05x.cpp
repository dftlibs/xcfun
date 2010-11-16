#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M05 exchange functional. to be used with HF exchange factor of .28 

template<class num>
static num ENERGY_FUNCTION(XC_M05X)(const densvars<num> &d)
{
   using m0xy_metagga_xc_internal::fw;

   const parameter param_a[12] =
     {  1.000000e+00,  8.151000e-02, -4.395600e-01, -3.224220e+00,  2.018190e+00,
        8.794310e+00, -2.950000e-03,  9.820290e+00, -4.823510e+00, -4.817574e+01,
        3.648020e+00,  3.402248e+01 };

   return (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua)
           + pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) );
}

NEW_TMGGA_FUNCTIONAL(XC_M05X);
SHORT_DESCRIPTION(XC_M05X) = "M05 exchange";
LONG_DESCRIPTION(XC_M05X) =
             "M05 Meta-Hybrid Exchange Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n";
TEST_VARS(XC_M05X) = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
TEST_MODE(XC_M05X) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_M05X) = 1;
TEST_THRESHOLD(XC_M05X) = 3e-5;
TEST_IN(XC_M05X) = {1., .8, 1., 1., 1., 0.165,   0.1050};
TEST_OUT(XC_M05X) =    { -1.57876583, -2.12127045, 
      -2.11264351, -0.00315462, 
      0.00000000,  -0.00444560,  
      3.45640232,   4.4349756 };

