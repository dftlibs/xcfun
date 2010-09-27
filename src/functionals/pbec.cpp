#include "functional.h"
#include "constants.h"
#include "pw92eps.h"
#include "vwn.h"


template<class num>
static num A(const num &eps, const num &u3)
{
  using xc_constants::param_beta_gamma;
  using xc_constants::param_gamma;
  return param_beta_gamma/expm1(-eps/(param_gamma*u3));
}

template<class num>
static num H(const num &d2, const num &eps, const num &u3)
{
  num d2A = d2*A(eps,u3);
  using xc_constants::param_gamma;
  using xc_constants::param_beta_gamma;
  return param_gamma*u3*
    log(1+param_beta_gamma*d2*(1 + d2A)/(1+d2A*(1+d2A)));
}

// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template<class num>
static num phi(const densvars<num> &d)
{
  return pow(2.0,-1.0/3.0)*d.n_m13*d.n_m13*(sqrt(d.a_43)+sqrt(d.b_43));
}

template<class num>
static num energy(const densvars<num> &d)
{
  num eps = pw92eps::pw92eps(d);
  num u = phi(d);
  // Avoiding the square root of d.gnn here
  num d2 = pow(1.0/12*pow(3,5.0/6.0)/pow(M_PI,-1.0/6),2)*
    d.gnn/(u*u*pow(d.n,7.0/3.0));
  return d.n*(eps + H(d2,eps,pow3(u)));
}

template<class num>
static num energy_vwn_pbec(const densvars<num> &d)
{
  num eps = vwn::vwn5_eps(d);
  num u = phi(d);
  // Avoiding the square root of d.gnn here
  num d2 = pow(1.0/12*pow(3,5.0/6.0)/pow(M_PI,-1.0/6),2)*
    d.gnn/(u*u*pow(d.n,7.0/3.0));
  return d.n*(eps + H(d2,eps,pow3(u)));
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
  const double d[5] = 
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
  const double ref[21] =
    {-0.184442072405E+01,   

     -0.814334534280E-01,
     -0.820182123795E-01,  
     0.510839298939E-06,
     0.102167859788E-05,
     0.510839298939E-06,  



     -0.124297349784E-02,
     -0.183505806584E-02,
     0.134850158624E-07, 
     0.269700317248E-07,  
     0.134850158624E-07,
     -0.125767116982E-02,
     0.136189478240E-07,
     0.272378956480E-07,    
     0.136189478240E-07,
     -0.216571369852E-12,
     -0.433142739704E-12,
     -0.216571369852E-12,  
     -0.866285479407E-12,
     -0.433142739704E-12,
     -0.216571369852E-12, 
    };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}


void setup_vwn_pbec(functional &f)
{
  f.describe(XC_VWN_PBEC, XC_GGA,
	     "PBE correlation functional using VWN LDA correlation.",
	     "PBE correlation functional with VWN LDA correlation.\n"
	     "J.P. Perdew, K. Burke, and M. Ernzerhof, Generalized\n"
	     "gradient approximation made simple, "
	     "Phys. Rev. Lett. 77 (1996) 3865-3868\n"
	     "This version of PBEc used VWN instead of PW92 as the LDA\n"
	     "correlation energy.\n"
	     "Implemented by Ulf Ekstrom\n");
  SET_GGA_ENERGY_FUNCTION(f,energy_vwn_pbec);
}
