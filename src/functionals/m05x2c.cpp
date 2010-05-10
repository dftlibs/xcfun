#include "functional.h"
#include "m0xy_fun.h"

// M05 correlation functional

template<class num>
static num energy (const densvars<num> &d)
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::m05_c_anti;
   using m0xy_metagga_xc_internal::ueg_c_anti;
   using m0xy_metagga_xc_internal::m05_c_para;
   using m0xy_metagga_xc_internal::ueg_c_para;

   // parameters for anti-parallel spin contributions
   static const parameter param_c_anti[5] =
       {  1.000000e+00,  1.092970e+00, -3.791710e+00,  2.828100e+00, -1.058909e+01 };

     // parameters for parallel spin contributions
   static const parameter param_c_para[5] =
       {  1.000000e+00, -3.054300e+00,  7.618540e+00,  1.476650e+00, -1.192365e+01 };

   num chi_a = chi(d.a, d.gaa);
   num chi_b = chi(d.b, d.gbb);
   num zet_a = zet(d.a, d.taua);
   num zet_b = zet(d.b, d.taub);

   num Ec_ab = ueg_c_anti(d)   * m05_c_anti(param_c_anti,chi_a,chi_b);
   num Ec_aa = ueg_c_para(d.a) * m05_c_para(param_c_para,chi_a,zet_a);
   num Ec_bb = ueg_c_para(d.b) * m05_c_para(param_c_para,chi_b,zet_b);

   return Ec_ab + Ec_aa + Ec_bb;
}

void setup_m05x2c(functional &f)
{
  f.describe("m05x2c",XC_MGGA,
	     "M05-2X Correlation",
             "M05-2X Meta-Hybrid Correlation Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  static const double d[] = 
    {1., .8, 1., 1., 1., .33, .21};
  static const double ref[] =
    { -0.06717000, -0.14727520,  0.04240607,  0.02498949,  0.03125835,  0.00000000, -0.07317847, -0.16011489 };
  f.add_test(XC_VARS_AB,1,d,ref,1e-5);
}
