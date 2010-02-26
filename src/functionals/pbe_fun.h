#include "pw92_fun.h"
#include "pw9xx_fun.h"

namespace PBEx_internal
{
   using pw91_like_x_internal::prefactor;

   // enhancement factor F(S), common to PBEx and REVPBEx
  template<class num>
  static num pbex_enhancement(const double &R, 
			      const num &rho,
			      const num &grad)
   {
      using pw91_like_x_internal::S;

      // double mu = 0.066725*M_PI*M_PI/3.0;
      // ulfek: mu from Daresbury implementation
      double mu = 0.2195149727645171;
      num st = S(rho,grad);
      num t1 = 1 + mu*st*st/R;
   
      return 1 + R - R/t1; 
   }

   // original PBE exchange functional

   template<class num>
   static num pbe_exchange(const densvars<num> &d)
   {
     const double R = 0.804;

     return (  prefactor(d.a,d.gaa)*pbex_enhancement(R,d.a,d.gaa) 
             + prefactor(d.b,d.gbb)*pbex_enhancement(R,d.b,d.gbb) );
   }

   // REVPBE exchange functional

   template<class num>
   static num revpbe_exchange(const densvars<num> &d)
   {
     const double R = 1.245;

     return (  prefactor(d.a,d.gaa)*pbex_enhancement(R,d.a,d.gaa) 
             + prefactor(d.b,d.gbb)*pbex_enhancement(R,d.b,d.gbb) );
   }
}


// PBEc implementation from MOLPRO manual, by ulfek.
namespace PBEc_internal
{
  //ulfek: I cannot find these constants directly in the PBE
  //paper
  //  const double iota = 0.0716;
  //MOLPRO value: const double kappa = 0.004235;
  //const double kappa = 0.00423488752945734;
  //const double nu = 16*pow(3*M_PI*M_PI,1.0/3)/M_PI;
  //const double lambda = nu*kappa;
  const double param_gamma = (1-log(2.0))/(M_PI*M_PI);
  //ulfek: beta comes from the slowly varying limit, see PBE paper
  const double param_beta_pbe_paper = 0.066725;
  const double param_beta_accurate  = 0.06672455060314922;
  //  const double param_omega = 0.046644;
  const double param_beta_gamma = param_beta_accurate/param_gamma;
  

  template<class num>
  static num u(const num &zeta)
  {
    // TODO: check if |zeta| ~= 1 
    return 0.5*(pow(1+zeta,2.0/3)+(pow(1-zeta,2.0/3)));
  }
  
  template<class num>
  static num A(const num &eps, const num &zeta)
  {
    return param_beta_gamma/(exp(-eps/(param_gamma*pow(u(zeta),3))) - 1);
  }

  template<class num>
  static num H(const num &d, const num &eps, const num &zeta)
  {
  num d2 = d*d;
  num d2A = d2*A(eps,zeta);
  return param_gamma*pow(u(zeta),3)*
    log(1+param_beta_gamma*d2*(1 + d2A)/(1+d2A*(1+d2A)));
  }

  template<class num>
  static num pbe_correlation(const densvars<num> &d)
  {
    const num &sigma = d.gnn;
    num dd = 1.0/12*(sqrt(sigma)*pow(3,5.0/6))/
      (u(d.zeta)*pow(M_PI,-1.0/6)*pow(d.n,7.0/6));
    num eps = pw92_eps_pbe(d);
    return d.n*(eps + H(dd,eps,d.zeta));
  }
}

template<>
struct functional<XC_PBE_CORRELATION>
{
  static const char *get_name(void) { return "pbec"; }
  static const char *get_reference(void) 
  { 
    return "PBE correlation functional.\n"
      "J.P. Perdew, K. Burke, and M. Ernzerhof, Generalized\n"
      "gradient approximation made simple, "
      "Phys. Rev. Lett. 77 (1996) 3865-3868";
  }
  enum { supported_modes = XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return PBEc_internal::pbe_correlation(d);
  }
  static int test(void) 
  { 
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_pbe.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      {-0.184442072405E+01,    -0.814334534280E-01,
       -0.820182123795E-01,     0.510839298939E-06,
        0.510839298939E-06,     0.102167859788E-05,
       -0.124297349784E-02,    -0.183505806584E-02,
        0.134850158624E-07,     0.134850158624E-07,
        0.269700317248E-07,    -0.125767116982E-02,
        0.136189478240E-07,     0.136189478240E-07,
        0.272378956480E-07,    -0.216571369852E-12,
       -0.216571369852E-12,    -0.433142739704E-12,
       -0.216571369852E-12,    -0.433142739704E-12,
       -0.866285479407E-12 };
    return standard_abgga_test<XC_PBE_CORRELATION>(d,ref,1e-11);
  }
};

template<>
struct functional<XC_PBE_EXCHANGE>
{
  static const char *get_name(void) { return "pbex"; }
  static const char *get_reference(void) 
  { 
    return "PBE exchange functional.\n"
      "J.P. Perdew, K. Burke, and M. Ernzerhof, Generalized\n"
      "gradient approximation made simple, "
      "Phys. Rev. Lett. 77 (1996) 3865-3868";    
  }
  enum { supported_modes = XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return PBEx_internal::pbe_exchange(d);
  }
  static int test(void) 
  { 
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_pbe.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      { -0.276589791995E+03,
	-0.382556082420E+01,
	-0.378108116179E+01,
	-0.174145337536E-04,
	-0.175120610339E-04,
	0.000000000000E+00,
	-0.429564214817E-01,
	0.000000000000E+00,
	0.185237729809E-06,
	0.000000000000E+00,
	0.000000000000E+00,
	-0.424802511645E-01,
	0.000000000000E+00,
	0.161839553501E-06,
	0.000000000000E+00,
	0.740514207206E-11,
	0.000000000000E+00,
	0.000000000000E+00,
	0.786563034093E-11,
	0.000000000000E+00,
	0.000000000000E+00
      };
    return standard_abgga_test<XC_PBE_EXCHANGE>(d,ref,1e-11);
  }
};

/*
template<>
struct functional<XC_REVPBE_EXCHANGE_UNTESTED>
{
  static const char *get_name(void) { return "revpbex"; }
  static const char *get_reference(void) { return "blah blah"; }
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return PBEx_internal::revpbe_exchange(d);
  }
  static int test(void) { return -1; }
};
*/

//TODO: deprecated
using PBEc_internal::pbe_correlation;
using PBEx_internal::pbe_exchange;
using PBEx_internal::revpbe_exchange;

