#include "pw9xx_fun.h"

namespace PW91xk_internal
{
   using pw91_like_x_internal::prefactor;
   using pw91_like_x_internal::pw91k_prefactor;

   // enhancement factor F(S), common to PW91x and PW91k 
   template<class num>
   static num pw91xk_enhancement(const double param_AB[6],
                                 const num &rho,
                                 const num &grad)
   {
      using pw91_like_x_internal::S;

      num st = S(rho,grad);

      num t1 = 1 + param_AB[0]*st*asinh(param_AB[1]*st); 
      num t2 = st*st*(param_AB[2] - param_AB[3]*exp(-param_AB[4]*st*st)); 

      num numerator   = t1 + t2; 
      num denominator = t1 + param_AB[5]*pow(st,4);

      return numerator/denominator;
   }

   // PW91 exchange functional

   template<class num>
   static num pw91_exchange(const num &rhoa,
                            const num &rhob,
                            const num &grada,
                            const num &gradb)
   {
     const static double param_AB[6] = 
       { 0.19645, 7.7956, 0.2743, 0.1508, 100.0, 0.004};

     return (  prefactor(rhoa,grada)*pw91xk_enhancement(param_AB,rhoa,grada)
             + prefactor(rhob,gradb)*pw91xk_enhancement(param_AB,rhob,gradb) );
   }

   // PW91 kinetic energy functional
   // reduces to TF kinetic energy functional for zero gradient

   template<class num>
   static num pw91_kinetic(const num &rhoa,
                           const num &rhob,
                           const num &grada,
                           const num &gradb)
   {
     const static double param_AB[6] = 
       { 0.093907, 76.320, 0.26608, 0.0809615, 100.0, 0.57767e-4};

     return (  pw91k_prefactor(rhoa)*pw91xk_enhancement(param_AB,rhoa,grada)
             + pw91k_prefactor(rhob)*pw91xk_enhancement(param_AB,rhob,gradb) );
   }
}

using PW91xk_internal::pw91_exchange;
using PW91xk_internal::pw91_kinetic;

