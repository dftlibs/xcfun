#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M06 exchange functional. to be used with HF exchange factor of .27 

template<class num>
static num energy (const densvars<num> &d)
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using m0xy_metagga_xc_internal::alpha_x;

   static const parameter param_a[12] = 
       {  5.877943e-01, -1.371776e-01,  2.682367e-01, -2.515898e+00, -2.978892e+00, 
          8.710679e+00,  1.688195e+01, -4.489724e+00, -3.299983e+01, -1.449050e+01,
          2.043747e+01,  1.256504e+01 }; 
   static const parameter param_d[6] = 
       {  1.422057e-01,  7.370319e-04, -1.601373e-02,  0.000000e+00,  0.000000e+00, 
          0.000000e+00 }; 

   num chia2 = pow(chi(d.a, d.gaa), 2);
   num chib2 = pow(chi(d.b, d.gbb), 2);

   return (  (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua) 
              + lsda_x(d.a)*h(param_d, alpha_x, chia2, zet(d.a, d.taua)))
            +
             (  pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) 
              + lsda_x(d.b)*h(param_d, alpha_x, chib2, zet(d.b, d.taub)))
          );
}

void setup_m06x(functional &f)
{
  f.describe("m06x",XC_MGGA,
	     "M06 exchange",
             "M06 Meta-Hybrid Exchange Functional\n"
             "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  static const double d[] =
    {1., .8, 1., 1., 1., .33, .21};
  static const double ref[] =
    { -0.93035619, -1.57676773, -1.59326350, -0.00081086, -0.00129163,  0.00000000,  1.62369085,  2.06708653 };
  f.add_test(XC_VARS_AB,1,d,ref,3e-5);
}
