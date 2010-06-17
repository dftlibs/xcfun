#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M06-2X exchange functional. to be used with HF exchange factor of .54
// because the param_d[5] array is all zero, it is not included here, and
// therefore the lsda_x() * h() terms (see e.g M06x) drop out, as h()=0

template<class num>
static num energy (const densvars<num> &d)
{
   using m0xy_metagga_xc_internal::fw;

   static const parameter param_a[12] =
     {  4.600000e-01, -2.206052e-01, -9.431788e-02,  2.164494e+00, -2.556466e+00,
       -1.422133e+01,  1.555044e+01,  3.598078e+01, -2.722754e+01, -3.924093e+01,
        1.522808e+01,  1.522227e+01 };

     return (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua)
             + pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) );
}

void setup_m06x2x(functional &f)
{
  f.describe(XC_M06X2X, XC_MGGA,
	     "M06-2X exchange",
             "M06-2X Meta-Hybrid Exchange Functional\n"
             "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  static const double d[] =
    {1., .8, 1., 1., 1., .33, .21};
  static const double ref[] =
    { -0.63803890, -0.81863653, -0.81208750, -0.00127795, -0.00179117,  0.00000000,  0.62610498,  0.80404158 }; 
  f.add_test(XC_VARS_AB,1,d,ref,3e-5);
}
