#include "functional.hpp"
#include "constants.hpp"
#include "pw92eps.hpp"


// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template<class num>
static num phi(const densvars<num> &d)
{
  return pow(2.0,-1.0/3.0)*d.n_m13*d.n_m13*(sqrt(d.a_43)+sqrt(d.b_43));
}

template<class num>
static num energy(const densvars<num> &d)
{
  using xc_constants::param_gamma;
  const parameter beta = 0.079030523241;
  num bg = beta/param_gamma;
  num eps = pw92eps::pw92eps(d);
  num u = phi(d);
  num u3 = pow3(u);
  // d2 is t^2
  num d2 = pow(1.0/12*pow(3,5.0/6.0)/pow(M_PI,-1.0/6),2.0)*
    d.gnn/(u*u*pow(d.n,7.0/3.0));
  num A = bg/expm1(-eps/(param_gamma*u3));
  num d2A = d2*A;
  num H = param_gamma*u3*
    log(1+bg*d2*(1 + d2A)/(1+d2A*(1+d2A)));
  return d.n*(eps + H);
}


FUNCTIONAL(XC_APBEC) = {
  "APBE correlation functional.",
  "APBE correlation functional\n"
  "L.A. Constantin, E. Fabiano, S. Laricchia, F. Della Sala,\n"
  "      Phys. Rev. Lett. 106, 186406 (2011).\n"
  "Implemented by Eduardo Fabiano\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(energy)
};
