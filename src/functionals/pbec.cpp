#include "functional.h"
#include "pw92eps.h"

template<class num>
static num u(const num &zeta)
{
  // TODO: check if |zeta| ~= 1 
  return 0.5*(pow(1+zeta,2.0/3.0)+(pow(1-zeta,2.0/3.0)));
}

template<class num>
static num A(const num &eps, const num &u3)
{
  using xc_constants::param_beta_gamma;
  using xc_constants::param_gamma;
  return param_beta_gamma/(exp(-eps/(param_gamma*u3)) - 1);
}

template<class num>
static num H(const num &d, const num &eps, const num &u3)
{
  num d2 = d*d;
  num d2A = d2*A(eps,u3);
  using xc_constants::param_gamma;
  using xc_constants::param_beta_gamma;
  return param_gamma*u3*
    log(1+param_beta_gamma*d2*(1 + d2A)/(1+d2A*(1+d2A)));
}

template<class num>
static num energy(const densvars<num> &d)
{
  num u1 = preexpand(u,d.zeta);
  num dd = 1.0/12*(sqrt(d.gnn)*pow(3,5.0/6.0))/
    (u1*pow(M_PI,-1.0/6)*pow(d.n,7.0/6.0));
  num eps = pw92eps::eps_pbe(d);
  return d.n*(eps + H(dd,eps,pow3(u1)));
}

void setup_pbec(functional &f)
{
  f.describe(XC_PBEC, XC_GGA,
	     "PBE correlation functional",
	     "PBE correlation functional.\n"
	     "J.P. Perdew, K. Burke, and M. Ernzerhof, Generalized\n"
	     "gradient approximation made simple, "
	     "Phys. Rev. Lett. 77 (1996) 3865-3868\n"
	     "Implemented by Ulf Ekstrom\n");
  SET_GGA_ENERGY_FUNCTION(f,energy);
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
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}
