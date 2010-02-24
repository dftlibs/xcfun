#include "m0xy_fun.h"


// Exchange functionals for the different M05 flavors 

namespace M05_x_internal
{
   using m0xy_metagga_xc_internal::fw;
   using PBEx_internal::pbe_exchange;

   // M05 exchange functional. to be used with HF exchange factor of .28 

   template<class num>
   static num m05_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  1.000000e+00,  8.151000e-02, -4.395600e-01, -3.224220e+00,  2.018190e+00, 
          8.794310e+00, -2.950000e-03,  9.820290e+00, -4.823510e+00, -4.817574e+01,
          3.648020e+00,  3.402248e+01 }; 

     return (  pbe_exchange(d.a,d.gaa)*fw(param_a, d.a, d.taua) 
             + pbe_exchange(d.b,d.gbb)*fw(param_a, d.b, d.taub) );
   }

   // M05-2X exchange functional. to be used with HF exchange factor of .56 

   template<class num>
   static num m05_2x_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  1.000000e+00, -5.683300e-01, -1.300570e+00,  5.500700e+00,  9.064020e+00, 
         -3.221075e+01, -2.373298e+01,  7.022996e+01,  2.988614e+01, -6.025778e+01,
         -1.322205e+01,  1.523694e+01 }; 

     return (  pbe_exchange(d.a,d.gaa)*fw(param_a, d.a, d.taua) 
             + pbe_exchange(d.b,d.gbb)*fw(param_a, d.b, d.taub) );
   }

}

using M05_x_internal::m05_exchange;
using M05_x_internal::m05_2x_exchange;


