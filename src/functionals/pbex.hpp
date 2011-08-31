#ifndef PBEX_H
#define PBEX_H
#include "pw9xx.hpp"

namespace pbex
{
  // values of R for PBEx and REVPBEx

  const parameter R_pbe = 0.804;
  const parameter R_revpbe = 1.245;

  // enhancement factor F(S), common to PBEx and REVPBEx
  template<class num>
  static num enhancement(const parameter &R, 
			 const num &rho,
			 const num &grad)
  {
    using pw91_like_x_internal::S2;
    
    // parameter mu = 0.066725*M_PI*M_PI/3.0;
    // ulfek: mu from Daresbury implementation
    const parameter mu = 0.2195149727645171;
    num st2 = S2(rho,grad);
    num t1 = 1 + mu*st2/R;   
    return 1 + R - R/t1; // Intel <= 11.1 miscompiles(?) this line with -fast
  }

  template<class num>
  static num enhancement_RPBE(const num &rho,
			      const num &grad)
  {
    using pw91_like_x_internal::S2;
    const parameter mu = 0.2195149727645171;
    return 1 - R_pbe*expm1((-mu/R_pbe)*S2(rho,grad));
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
