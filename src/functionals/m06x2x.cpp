#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M06-2X exchange functional. to be used with HF exchange factor of .54
// because the param_d[5] array is all zero, it is not included here, and
// therefore the lsda_x() * h() terms (see e.g M06x) drop out, as h()=0

template<class num>
static num ENERGY_FUNCTION(XC_M06X2X)(const densvars<num> &d)
{
   using m0xy_metagga_xc_internal::fw;

   const parameter param_a[12] =
     {  4.600000e-01, -2.206052e-01, -9.431788e-02,  2.164494e+00, -2.556466e+00,
       -1.422133e+01,  1.555044e+01,  3.598078e+01, -2.722754e+01, -3.924093e+01,
        1.522808e+01,  1.522227e+01 };

     return (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua)
             + pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) );
}

NEW_TMGGA_FUNCTIONAL(XC_M06X2X);
SHORT_DESCRIPTION(XC_M06X2X) = "M06-2X exchange";
LONG_DESCRIPTION(XC_M06X2X) =
             "M06-2X Meta-Hybrid Exchange Functional\n"
             "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
             "Implemented by Andre Gomes\n";
TEST_VARS(XC_M06X2X) = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
TEST_MODE(XC_M06X2X) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_M06X2X) = 1;
TEST_THRESHOLD(XC_M06X2X) = 3e-5;
TEST_IN(XC_M06X2X) = {1., .8, 1., 1., 1., 0.165,   0.1050};
TEST_OUT(XC_M06X2X) =
    { -0.63803890, -0.81863653, -0.81208750, -0.00127795,  0.00000000,  -0.00179117, 1.25220996,   1.60808316 }; 

