#include "functional.hpp"
#include "pw9xx.hpp"

template<class num>
static num apbe_enhancement(const num &na, const num &gaa)
{
  using pw91_like_x_internal::S2;
  const parameter mu = 0.26;
  const parameter kappa = 0.804;
  num st2 = S2(na,gaa);
  num t1 = 1 + mu*st2/kappa;
  return 1 + kappa - kappa/t1;
}


template<class num>
static num energy_apbex(const num &na, const num &gaa)
{
  const parameter c = pow(81/(4*M_PI),1.0/3.0)/2;
  num na43 = pow(na,4.0/3.0);
  num lda = -c*na43;
  num apbex = lda * apbe_enhancement(na,gaa);
  return apbex;
}

template<class num>
static num energy(const densvars<num> &d)
{
  return energy_apbex(d.a,d.gaa) + energy_apbex(d.b,d.gbb);
}


FUNCTIONAL(XC_APBEX) = {
  "APBE Exchange Functional",
  "APBE Exchange Functional\n"
  "L.A. Constantin, E. Fabiano, S. Laricchia, F. Della Sala,\n"
  "      Phys. Rev. Lett. 106, 186406 (2011).\n"
  "Implemented by Eduardo Fabiano\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energy)
};

