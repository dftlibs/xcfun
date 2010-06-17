#include "functional.h"
#include "pw9xx.h"

// PW91 kinetic energy functional
// reduces to TF kinetic energy functional for zero gradient

template<class num>
static num energy(const densvars<num> &d)
{
  static const parameter param_AB[6] = 
    { 0.093907, 76.320, 0.26608, 0.0809615, 100.0, 0.57767e-4};
  using pw91_like_x_internal::pw91k_prefactor;
  using pw91_like_x_internal::pw91xk_enhancement;
  return pw91k_prefactor(d.a)*pw91xk_enhancement(param_AB,d.a,d.gaa)
    + pw91k_prefactor(d.b)*pw91xk_enhancement(param_AB,d.b,d.gbb);
}

void setup_pw91k(functional &f)
{
  f.describe(XC_PW91K, XC_LDA,
	     "PW91 GGA Kinetic Energy Functional",
	     "PW91 GGA Kinetic Energy Functional\n"
	     "A. Lembarki, H. Chermette, Phys. Rev. A 50, 5328 (1994)\n"
	     "Implemented by Andre Gomes.\n"); 
  SET_LDA_ENERGY_FUNCTION(f,energy);
}
