#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M05-2X exchange functional. to be used with HF exchange factor of .56 

template<class num>
static num energy (const densvars<num> &d)
{
   using m0xy_metagga_xc_internal::fw;

   static const parameter param_a[12] =
     {  1.000000e+00, -5.683300e-01, -1.300570e+00,  5.500700e+00,  9.064020e+00,
       -3.221075e+01, -2.373298e+01,  7.022996e+01,  2.988614e+01, -6.025778e+01,
       -1.322205e+01,  1.523694e+01 };

   return (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua)
           + pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) );
}

void setup_m05x2x(functional &f)
{
  f.describe(XC_M05X2X,XC_MGGA,
	     "M05-2X exchange",
             "M05-2X Meta-Hybrid Exchange Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  static const double d[] =
    {1., .8, 1., 1., 1., .33, .21};
  static const double ref[] =
    { -1.38233309, -0.19638222, -0.08614105, -0.00289174, -0.00365982,  0.00000000, -1.59421158, -1.95293869 }; 
  f.add_test(XC_VARS_AB,1,d,ref,3e-5);
}
