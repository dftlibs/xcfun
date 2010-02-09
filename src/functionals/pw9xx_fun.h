#ifndef PW9xx_H
#define PW9xx_H

// common functions for exchange (and kinetic energy) functionals
//
// related to PW91x - like PW91x (and PW91k, which is a reparametrization 
// of PW91x to describe T_s), PBEx, PBEREVx and so on 

namespace pw91_like_x_internal
{

  const char *pw91x_reference = 
    "Perdew-Wang 1991 GGA Exchange Functional\n"
      "J. P. Perdew, J. A. Chevary, S. H. Vosko, "
      "K. A. Jackson, M. R. Pederson, and C. Fiolhais, "
      "Phys. Rev. B 46, 6671 (1992)\n"; 

  const char *pw91k_reference = 
    "PW91 GGA Kinetic Energy Functional\n"
      "A. Lembarki, H. Chermette, "
      "Phys. Rev. A 50, 5328 (1994)\n"; 

  const char *pbex_reference = 
    "PBE Exchange Functional\n"
      "J. P. Perdew, K. Burke, and M. Ernzerhof, "
      "Phys. Rev. Lett 77, 3865 (1996)\n"; 


  const char *pberevx_reference = 
    "Revised PBE Exchange Functional\n"
      "Y. Zhang and W., "
      "Phys. Rev. Lett 80, 890 (1998)\n";


// formulas for the auxiliary quantities below can be found, for
// instance, at 
//
//    http://www.molpro.net/info/current/doc/manual/node734.html#dftfun:PW91X
//

  template<class num>
  static num chi(const num &rho, const num &grad)
  {
    return sqrt(grad)/pow(rho,4.0/3.0);
  }

  template<class num>
  static num S(const num &rho, const num &grad)
  {
    return chi(rho,grad)*pow(6.0,2.0/3.0)/(12*pow(M_PI,2.0/3.0)); 
  }

// prefactor multiples the enhancement factor F(S), which is then different
// for the different functionals

  template<class num>
  static num prefactor(const num &rho, const num &grad)
  {
// aspg: the 2^.333 factor here i can't see in the molpro formula, will have to 
// double-check this - but as it is it matches the results for the database
// e.g
//   http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_pw91.html
//   http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_pbe.html
// and for pbe, it is close to the dalton code..(trunncation of variables
// there make things deviate?)
// ue: More likely the parameters, Dalton has some strange choice of
// precision in its constants.
    return -0.75*pow(2.0,1.0/3.0)*pow(3*M_PI*M_PI,1.0/3.0)
      *pow(rho,4.0/3.0)/M_PI;
  }

// prefactor for the pw91k functional
  template<class num>
  static num pw91k_prefactor(const num &rho)
  {
    using xc_constants::CF;

    return CF*pow(2.0,2.0/3.0)*pow(rho,5.0/3.0);
  }

}

#endif
