#include "functional.h"
#include "m0xy_fun.h"

// M05 correlation functional

template<class num>
static num ENERGY_FUNCTION(XC_M05C)(const densvars<num> &d)
{
   using pw91_like_x_internal::chi2;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::m05_c_anti;
   using m0xy_metagga_xc_internal::ueg_c_anti;
   using m0xy_metagga_xc_internal::m05_c_para;
   using m0xy_metagga_xc_internal::ueg_c_para;

   // parameters for anti-parallel spin contributions
   const parameter param_c_anti[5] =
     {  1.000000e+00,  3.785690e+00, -1.415261e+01, -7.465890e+00,  1.794491e+01 };

   // parameters for parallel spin contributions
   const parameter param_c_para[5] =
     {  1.000000e+00,  3.773440e+00, -2.604463e+01,  3.069913e+01, -9.226950e+00 };

   num chi_a2 = chi2(d.a, d.gaa);
   num chi_b2 = chi2(d.b, d.gbb);
   num zet_a = zet(d.a, d.taua);
   num zet_b = zet(d.b, d.taub);
   num Dsigma_a = m0xy_metagga_xc_internal::Dsigma(d.a,d.gaa,d.taua);
   num Dsigma_b = m0xy_metagga_xc_internal::Dsigma(d.b,d.gbb,d.taub);

   num Ec_ab = ueg_c_anti(d)   * m05_c_anti(param_c_anti,chi_a2,chi_b2);
   num Ec_aa = ueg_c_para(d.a) * m05_c_para(param_c_para,chi_a2,zet_a, Dsigma_a);
   num Ec_bb = ueg_c_para(d.b) * m05_c_para(param_c_para,chi_b2,zet_b, Dsigma_b);

   return Ec_ab + Ec_aa + Ec_bb;
}

NEW_TMGGA_FUNCTIONAL(XC_M05C);
SHORT_DESCRIPTION(XC_M05C) = "M05 Correlation";
LONG_DESCRIPTION(XC_M05C) =
	     "M05 Meta-Hybrid Correlation Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n";
NO_TEST(XC_M05C);
  /*  const double d[] = 
    {1., .8, 1., 1., 1., .33, .21};
  const double ref[] =
    { -0.06599246, -0.15418438,  0.02729798,  0.03090769,  0.04788146,  0.00000000, -0.09757618, -0.23742358 };
  f.add_test(XC_VARS_AB,1,d,ref,2e-5); */

