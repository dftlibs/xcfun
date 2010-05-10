#ifndef PBEX_H
#define PBEX_H
#include "pw9xx.h"

namespace pbex
{
  // values of R for PBEx and REVPBEx

  static const parameter R_pbe = 0.804;
  static const parameter R_revpbe = 1.245;

  // enhancement factor F(S), common to PBEx and REVPBEx
  template<class num>
  static num enhancement(const parameter &R, 
			 const num &rho,
			 const num &grad)
  {
    using pw91_like_x_internal::S;
    
    // parameter mu = 0.066725*M_PI*M_PI/3.0;
    // ulfek: mu from Daresbury implementation
    static const parameter mu = 0.2195149727645171;
    num st = S(rho,grad);
    num t1 = 1 + mu*st*st/R;   
    return 1 + R - R/t1; // Intel <= 11.1 miscompiles(?) this line with -fast
  }

  template<class num>
  static num energy_pbe_ab(const parameter &R,
			   const num &rho,
			   const num &grad)
  {
    using pw91_like_x_internal::prefactor;
    return prefactor(rho)*enhancement(R,rho,grad);
  }
}

#endif
