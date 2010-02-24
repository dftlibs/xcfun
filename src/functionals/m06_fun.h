#include "m0xy_fun.h"

namespace M06_x_internal
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using PBEx_internal::pbe_exchange;

   // M06-L exchange functional. Does not use HF exchange. 

   template<class num>
   static num m06_l_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  3.987756e-01,  2.548219e-01,  2.923994e-01, -2.103655e+00, -6.302147e+00, 
          1.097615e+01,  3.097273e+01, -2.318489e+01, -5.673480e+01,  2.160364e+01,
          3.421814e+01, -9.049762e+00 }; 
     const static double param_d[5] = 
       {  6.012244e-01,  4.748822e-03, -8.635108e-04, -9.308062e-06,  4.482811e-05,
          0.000000e+00 }; 

     return (  (  pbe_exchange(d.a, d.gaa)*fw(param_a, d.a, d.taua) 
                + lsda_x(d.a)*h(param_d, chi(d.a, d.gaa), zet(d.a, d.taua)))
              +
               (  pbe_exchange(d.b, d.gbb)*fw(param_a, d.b, d.taub) 
                + lsda_x(d.b)*h(param_d, chi(d.b, d.gbb), zet(d.b, d.taub)))
            );
   } 

   // M06-HF exchange functional. Uses full HF exchange. 

   template<class num>
   static num m06_hf_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  1.179732e-01, -1.066708e+00, -1.462405e-01,  7.481848e+00,  3.776679e+00, 
         -4.436118e+01, -1.830962e+01,  1.003903e+02,  3.864360e+01, -9.806018e+01,
          2.557716e+01,  3.590404e+01 }; 
     const static double param_d[5] = 
       { -1.179732e-01, -2.500000e-03, -1.180056e-02,  0.000000e+00,  0.000000e+00,
          0.000000e+00 }; 

     return (  (  pbe_exchange(d.a, d.gaa)*fw(param_a, d.a, d.taua) 
                + lsda_x(d.a)*h(param_d, chi(d.a, d.gaa), zet(d.a, d.taua)))
              +
               (  pbe_exchange(d.b, d.gbb)*fw(param_a, d.b, d.taub) 
                + lsda_x(d.b)*h(param_d, chi(d.b, d.gbb), zet(d.b, d.taub)))
            );
   } 

   // M06 exchange functional. to be used with HF exchange factor of .27 

   template<class num>
   static num m06_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  5.877943e-01, -1.371776e-01,  2.682367e-01, -2.515898e+00, -2.978892e+00, 
          8.710679e+00,  1.688195e+01, -4.489724e+00, -3.299983e+01, -1.449050e+01,
          2.043747e+01,  1.256504e+01 }; 
     const static double param_d[5] = 
       {  1.422057e-01,  7.370319e-04, -1.601373e-02,  0.000000e+00,  0.000000e+00,
          0.000000e+00 }; 

     return (  (  pbe_exchange(d.a, d.gaa)*fw(param_a, d.a, d.taua) 
                + lsda_x(d.a)*h(param_d, chi(d.a, d.gaa), zet(d.a, d.taua)))
              +
               (  pbe_exchange(d.b, d.gbb)*fw(param_a, d.b, d.taub) 
                + lsda_x(d.b)*h(param_d, chi(d.b, d.gbb), zet(d.b, d.taub)))
            );
   }

   // M06-2X exchange functional. to be used with HF exchange factor of .54 
   // because the param_d[5] array is all zero, it is not included here, and
   // therefore the lsda_x() * h() terms (see above) drop out, as h()=0

   template<class num>
   static num m06_2x_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  4.600000e-01, -2.206052e-01, -9.431788e-02,  2.164494e+00, -2.556466e+00, 
         -1.422133e+01,  1.555044e+01,  3.598078e+01, -2.722754e+01, -3.924093e+01,
          1.522808e+01,  1.522227e+01 }; 

     return (  pbe_exchange(d.a, d.gaa)*fw(param_a, d.a, d.taua) 
             + pbe_exchange(d.b, d.gbb)*fw(param_a, d.b, d.taub) );
   } 

}

using M06_x_internal::m06_l_exchange;
using M06_x_internal::m06_hf_exchange;
using M06_x_internal::m06_exchange;
using M06_x_internal::m06_2x_exchange;

