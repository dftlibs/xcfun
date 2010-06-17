#include "functional.h"
#include "m0xy_fun.h"

// M06 correlation functional

template<class num>
static num energy (const densvars<num> &d)
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::m06_c_anti;
   using m0xy_metagga_xc_internal::ueg_c_anti;
   using m0xy_metagga_xc_internal::m06_c_para;
   using m0xy_metagga_xc_internal::ueg_c_para;

   // parameters for anti-parallel spin contributions
   static const parameter param_c_anti[5] =
     {  3.741539e+00,  2.187098e+02, -4.531252e+02,  2.936479e+02, -6.287470e+01 };
   static const parameter param_d_anti[6] =
     { -2.741539e+00, -6.720113e-01, -7.932688e-02,  1.918681e-03, -2.032902e-03, 
        0.000000e+00 };

   // parameters for parallel spin contributions
   static const parameter param_c_para[5] =
     {  5.094055e-01, -1.491085e+00,  1.723922e+01, -3.859018e+01,  2.845044e+01 };
   static const parameter param_d_para[6] =
     {  4.905945e-01, -1.437348e-01,  2.357824e-01,  1.871015e-03, -3.788963e-03,
        0.000000e+00 }; 

   num chi_a = chi(d.a, d.gaa);
   num chi_b = chi(d.b, d.gbb);
   num zet_a = zet(d.a, d.taua);
   num zet_b = zet(d.b, d.taub);

   num Ec_ab = ueg_c_anti(d) * m06_c_anti(param_c_anti,param_d_anti,chi_a,zet_a, 
                                                                    chi_b,zet_b);
   num Ec_aa = ueg_c_para(d.a) * m06_c_para(param_c_para,param_d_para,chi_a,zet_a);
   num Ec_bb = ueg_c_para(d.b) * m06_c_para(param_c_para,param_d_para,chi_b,zet_b);

   return Ec_ab + Ec_aa + Ec_bb;
}

void setup_m06c(functional &f)
{
  f.describe(XC_M06C, XC_MGGA,
	     "M06 Correlation",
             "M06 Meta-Hybrid Correlation Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  static const double d[] = 
    {1., .8, 1., 1., 1., .33, .21};
  static const double ref[] =
    { -1.57876583, -2.12127045, -2.11264351, -0.00315462, -0.00444560,  0.00000000,  1.72820116,  2.21748787 };
  f.add_test(XC_VARS_AB,1,d,ref,1e-5);
}

