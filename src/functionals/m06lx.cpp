#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M06-L exchange functional. does not use HF exchange 

template<class num>
static num energy (const densvars<num> &d)
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using m0xy_metagga_xc_internal::alpha_x;

     static const parameter param_a[12] =
       {  3.987756e-01,  2.548219e-01,  3.923994e-01, -2.103655e+00, -6.302147e+00,
          1.097615e+01,  3.097273e+01, -2.318489e+01, -5.673480e+01,  2.160364e+01,
          3.421814e+01, -9.049762e+00 };
     static const parameter param_d[6] =
       {  6.012244e-01,  4.748822e-03, -8.635108e-03, -9.308062e-06,  4.482811e-05,
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

void setup_m06lx(functional &f)
{
  f.describe("m06lx",XC_MGGA,
	     "M06-L exchange",
             "M06-L Meta-GGA Exchange Functional\n"
             "Y Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101 (2006)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);
  static const double d[] =
    {1., .8, 1., 1., 1., .33, .21};
  static const double ref[] =
    { -1.60059999, -1.85109109, -1.81344820, -0.00370516, -0.00508827,  0.00000000,  1.19430735,  1.52656126 }; 
  f.add_test(XC_VARS_AB,1,d,ref,1e-5);
}
