#ifndef M0xy_H
#define M0xy_H

#include "pw92eps.h"
#include "pw9xx.h"
#include <stdio.h>

// common functions for MO5 and M06 family of (hybrid) meta-gga functionals

namespace m0xy_metagga_xc_internal
{

  // additional information from the authors, and reference implementations 
  // in F77, can be found at
  //
  //    http://comp.chem.umn.edu/mfm/

  // values for alpha are taken from
  //   T. Van Voorhis and G. E. Scuseria, J. Chem. Phys. 129, 219901 (2008)
  //   DOI: 10.1063/1.3005348

  static const parameter alpha_x = 0.00186726;
  static const parameter alpha_c_parallel = 0.00515088;
  static const parameter alpha_c_antiparallel = 0.00304966;


  // scalefactorTFconst is here because in the m05/m06/m08 papers
  // C_F is not your usual thomas-fermi constant, (3/10)*(3 pi^2)^(2/3)
  // but rather (3/5)*(6 pi^2)^(2/3) so a factor of (2^(5/2))^(2/3) = 3.174802 
  // is needed to multiply C_F
 
  static const parameter scalefactorTFconst = 3.174802;

  // zet calculates one of the working variables
  // see see TCA (2008) 120:215 eq 3
  // rho is the density, tau the kinetic energy density
  // 
  // andre gomes, 2010-04-14:
  //
  // one here must be aware that the factor of two is *not* in the
  // actual formulae in the paper - eq 3 - but is due to the use of
  // a definition of tau that is twice the value of spin kinetic energy
  // densities, so that tau -> tau/2 and therefore 2tau/rho53 -> tau/rho53.
  //
  // for further information one should consult the documentation
  // of the minnesota functional module 1.2 (MFM1.2) and the actual MFM code,
  // which was used as benchmark here for energy and first derivatives.

  template<class num>
  static num zet(const num &rho, const num &tau)
  {
    using xc_constants::CF;

    return ((tau/pow(rho,5.0/3.0)) - CF*scalefactorTFconst);
  }

  // gamma calculates the quantity of TCA (2008) 120:215 eq 4, but we supply it with chi^2 as parameter

  template<class num>
  static num gamma(const parameter alpha, const num &chi2, const num &zet)
  {
    return 1 + alpha*(chi2 + zet); 
  }

  // h calculates the quantity of TCA (2008) 120:215 eq 5, but here we provide it with the square of the x_\sigma
  // variable, since that is used always in terms of x^2. the d[] parameter array should be 
  // of dimension 6..

  template<class num>
  static num h(const parameter d[], const parameter alpha, const num &chi2, const num &zet)
  {
    num gam1 = gamma(alpha,chi2,zet);

    num t1 =  d[0]/gam1;
    num t2 = (d[1]*chi2 + d[2]*zet)/(gam1 * gam1);
    num t3 = (d[3]*chi2*chi2 + d[4]*chi2*zet + d[5]*zet*zet)/(gam1 * gam1 * gam1);

    return t1+t2+t3;
  }


  // functions specific to the exchange functinals

  // fw is the spin kinetic enrgy density enhancement factor, see e.g.
  // TCA 2008 120:215 eq 8
  // here the number of parameters is hardcoded since all M05/M06 functionals
  // described at the time of coding (feb 2010) use 12 of them.

  template<class num>
  static num fw(const parameter a[12], const num &rho, const num &tau)
  {
    using pw91_like_x_internal::pw91k_prefactor;

    // andre gomes, 2010-04-14: 
    //
    // one here must be aware that the factor of two is *not* in the
    // actual formulae in the paper - eq 11 - but is due to the use of 
    // a definition of tau that is twice the value of spin kinetic energy
    // densities, so that tau -> tau/2 . 
    //
    // for further information one should consult the documentation 
    // of the minnesota functional module 1.2 (MFM1.2) and the actual MFM code,
    // which was used as benchmark here for energy and first derivatives.
 
    num tau_lsda = 2*pw91k_prefactor(rho);

    num  t = tau_lsda/tau;
    num  w = (t - 1)/(t + 1); 
    num fw = 0.0;

    for (int i = 0; i < 12; i++)
       fw += a[i] * pow(w,i);

    return fw;
  }


  // functions specific to the correlation functionals

  // Dsigma is the self-interaction correction factor, for parallel spins, used in M06 and M06-2X 
  // see for instance TCA 2008 120:215 eq 16
  // obviously, chi and zet are to have same spin: Dsigma(x_\alpha, z_\alpha) or Dsigma(x_\beta, z_\beta)
  // chi2 is chi*chi
  //

  template<class num>
  static num Dsigma(const num &chi2, const num &zet)
  {
    using xc_constants::CF;

    return (1.0 - 0.25*chi2/(zet + CF*scalefactorTFconst)); 
  }

  // g is an auxiliary function that can be used in connection to both
  // parallel and antiparallel spins (provided in the latter that the variable
  // passed on is not chi_\sigma^2 but (chi_a^2 + chi_b^2)). 
  //
  // for compactness, we also pass (chi_\sigma^2 or (chi_a^2 + chi_b^2)) multiplied
  // by the corresponding gamma factor, since the two are always multiplied together.
  //
  // see TCA 2008 120:215 eqs 13 and 15 for the actual formulas
  //
  // here the number of parameters is hardcoded since all M05/M06 functionals
  // described at the time of coding (feb 2010) use 5 of them.

  template<class num>
  static num g(const parameter param_c[5], const num &gamma_chi_squared)
  {
     num g = 0.0;
     num b   = gamma_chi_squared / (1.0 + gamma_chi_squared); 
  
     for (int i = 0; i < 5; i++)
       g += param_c[i] * pow(b,i);

    return g;
  }


  // m06_c_anti calculates the  (g + h) part of the integrand, for the M06y correlationg
  // functionals, for the antiparallel spins
  //
  // m06_c_para calculates the  (g + h)D part of the integrand, for the M06y correlationg
  // functionals, for parallel spins

  template<class num>
  static num m06_c_anti(const parameter param_c[5], 
                        const parameter param_d[5], 
                        const num &chi_a, const num &zet_a, 
                        const num &chi_b, const num &zet_b) 
  {
    const parameter gamma_c_anti = 0.0031;  // this is an "universal" constant for all M05/M06 functionals 

    const num zet_ab  = zet_a + zet_b;
    const num chi_ab2 = chi_a*chi_a + chi_b*chi_b;
    return g(param_c,gamma_c_anti*chi_ab2) + h(param_d,alpha_c_antiparallel,chi_ab2,zet_ab);
  }

  template<class num>
  static num m06_c_para(const parameter param_c[5], 
                        const parameter param_d[5], 
                        const num &chi, const num &zet) 
  {
    const parameter gamma_c_parallel = 0.06;  // this is an "universal" constant for all M05/M06 functionals 

    const num chi2 = chi*chi;

    return (g(param_c,gamma_c_parallel*chi2) + h(param_d,alpha_c_parallel,chi,zet))*Dsigma(chi2,zet);
  }


  // m05_c_anti calculates the  g part of the integrand, for the M05y correlationg
  // functionals, for the antiparallel spins
  //
  // m05_c_para calculates the  gD part of the integrand, for the M05y correlationg
  // functionals, for parallel spins

  template<class num>
  static num m05_c_anti(const parameter param_c[5], const num &chi_a, const num &chi_b) 
  {
    const parameter gamma_c_anti = 0.0031;  // this is an "universal" constant for all M05/M06 functionals 

    const num chi_ab2 = chi_a*chi_a + chi_b*chi_b;

    return g(param_c,gamma_c_anti*chi_ab2);
  }

  template<class num>
  static num m05_c_para(const parameter param_c[5], const num &chi, const num &zet) 
  {
    const parameter gamma_c_parallel = 0.06;  // this is an "universal" constant for all M05/M06 functionals 

    const num chi2 = chi*chi;

    return g(param_c,gamma_c_parallel*chi2)*Dsigma(chi2,zet);
  }


  // functions to calculate the \epsilon^{UEC}^{[\alpha\beta,\sigma\sigma]} factor in the
  // M05/M06 (and VSXC) functionals. The definitions of those can be found in 
  // 
  // Zhao, Schultz, Truhlar, J. Chem. Theory Comput. 2006, 2, 264   equations 7,8
  //  

  template<class num>
  static num ueg_c_anti(const densvars<num> &d)
  {
    using pw92eps::eps;
    return (eps(d)*d.n) - ueg_c_para(d.a) - ueg_c_para(d.b);
  }

  template<class num>
  static num ueg_c_para(const num &rho)
  {
    return pw92eps::eps_polarized_m0y(rho)*rho;
  }


}

// local spin density approximation for the exchange
template<class num>
static num lsda_x(const num &rho)
{
  return -(3.0/2.0)*pow(3.0/(4.0*M_PI), 1.0/3.0)*pow(rho, 4.0/3.0);
}


#endif
