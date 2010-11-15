#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M06 exchange functional. to be used with HF exchange factor of .27 

template<class num>
static num ENERGY_FUNCTION(XC_M06X)(const densvars<num> &d)
{
   using pw91_like_x_internal::chi2;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using m0xy_metagga_xc_internal::alpha_x;

   const parameter param_a[12] = 
       {  5.877943e-01, -1.371776e-01,  2.682367e-01, -2.515898e+00, -2.978892e+00, 
          8.710679e+00,  1.688195e+01, -4.489724e+00, -3.299983e+01, -1.449050e+01,
          2.043747e+01,  1.256504e+01 }; 
   const parameter param_d[6] = 
       {  1.422057e-01,  7.370319e-04, -1.601373e-02,  0.000000e+00,  0.000000e+00, 
          0.000000e+00 }; 

   num chia2 = chi2(d.a, d.gaa);
   num chib2 = chi2(d.b, d.gbb);

   return (  (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua) 
              + lsda_x(d.a)*h(param_d, alpha_x, chia2, zet(d.a, d.taua)))
            +
             (  pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) 
              + lsda_x(d.b)*h(param_d, alpha_x, chib2, zet(d.b, d.taub)))
          );
}

NEW_TMGGA_FUNCTIONAL(XC_M06X);
SHORT_DESCRIPTION(XC_M06X) = "M06 exchange";
LONG_DESCRIPTION(XC_M06X) =
             "M06 Meta-Hybrid Exchange Functional\n"
             "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
             "Implemented by Andre Gomes\n"
	     "Reference gradient from ADF\n";

  /*
  const double d[] =
    {1., .8, 1., 1., 1., .33, .21};
  const double ref[] =
    { -0.93035619, -1.57676773, -1.59326350, -0.00081086,   0.00000000, -0.00129163, 1.62369085,  2.06708653 };
  */
TEST_VARS(XC_M06X) = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
TEST_ORDER(XC_M06X) = 1;
TEST_THRESHOLD(XC_M06X) = 1e-7;
TEST_IN(XC_M06X) =  {0.153652558932587,
			     0.153652558932587,     
			     8.390981882024769E-002,
			     8.390981882024769E-002, 
			     8.390981882024769E-002,			     
			     6.826262722466429E-002,  
			     6.826262722466429E-002};
TEST_OUT(XC_M06X) = {-1.0989923683183833e-01, // Energy from xcfun
			       -0.519568215542419,
			       -0.519568215542419,
			       -1.757089749188437E-002,
			       0.000000000000000E+000,
			       -1.757089749188437E-002,
			       9.227728385689479E-002,
			       9.227728385689479E-002};

