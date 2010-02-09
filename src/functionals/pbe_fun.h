#include "pw92_fun.h"
#include "pw9xx_fun.h"

template<class num>
static num rpbe_exchange_dalton(const num &rhoa,
                         const num &rhob,
                         const num &grada,
                         const num &gradab,
                         const num &gradb)
{
  double t1 = pow(2.0,.3333333333333333);

   /* code */
  return 0.5*(-1.477117532764045*(0.804*(1.0-1.0/exp(.004492799480231091*pow(gradb,2.0)/pow(rhob,2.666666666666667)))+
        1.0)*t1*pow(rhob,1.333333333333333)-1.477117532764045*(0.804*
        (1.0-1.0/exp(.004492799480231091*pow(grada,
        2.0)/pow(rhoa,2.666666666666667)))+1.0)*t1*pow(rhoa,1.333333333333333));
}

template<class num>
static num pbe_exchange_dalton(const num &rhoa,
                            const num &rhob,
                            const num &grada2,
                            const num &gradb2)
{
  double t1 = pow(2.0,0.33333333333333);
  // TODO: test if rhoa or rhob are too small
  return  0.5*(-1.477117532764045*t1*pow(rhob,1.333333333333333)*
             (1.804-0.804/(0.00449279948023*gradb2/pow(rhob,2.666666666666667)+
                           1.0))
               -1.477117532764045*t1*pow(rhoa,1.333333333333333)*
             (1.804-0.804/(0.00449279948023*grada2/pow(rhoa,2.666666666666667)+
             1.0)));
}

namespace PBEx_internal
{
   using pw91_like_x_internal::prefactor;

   // enhancement factor F(S), common to PBEx and REVPBEx
   template<class num>
   static num pbex_enhancement(const num &R, 
                               const num &rho,
                               const num &grad)
   {
      using pw91_like_x_internal::S;

      num mu = 0.066725*M_PI*M_PI/3.0;
      num st = S(rho,grad);
      num t1 = 1.0 + (mu*st*st/R);
   
      // TODO: veryfy 1/t1 does not overflow/underflow

      return (1.0 + R - (R/t1)); 
   }

   // original PBE exchange functional

   template<class num>
   static num pbe_exchange(const num &rhoa,
			   const num &rhob,
			   const num &grada,
			   const num &gradb)
   {
     static num R = 0.804;

     return (  prefactor(rhoa,grada)*pbex_enhancement(R,rhoa,grada) 
             + prefactor(rhob,gradb)*pbex_enhancement(R,rhob,gradb) );
   }

   // REVPBE exchange functional

   template<class num>
   static num revpbe_exchange(const num &rhoa,
   	    	     	      const num &rhob,
			      const num &grada,
			      const num &gradb)
   {
     static num R = 1.245;

     return (  prefactor(rhoa,grada)*pbex_enhancement(R,rhoa,grada) 
             + prefactor(rhob,gradb)*pbex_enhancement(R,rhob,gradb) );
   }
}


// PBE correlation

template<class num>
static num pbe_correlation_unpolarized(const num &R,
				       const num &Z)
{
  num t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

  t1 = R;
  t2 = 1/pow(t1,0.33333333333333);
  t3 = 0.1325688999052*t2+1.0;
  t4 = log(16.0818243221511/(2.225569421150687*t2+0.18970043257476/
      pow(t1,0.66666666666667)+5.98255043577108/pow(t1,0.16666666666667)+
        0.80042863499936/sqrt(t1))+1.0);
  t5 = Z;
  t6 = 1/pow(t1,2.333333333333334);
  t7 = t5*t5;
  t8 = 1/pow(t1,4.666666666666667);
  t9 = exp(1.999923215200099*t3*t4)-1.0;
  t10 = 1/t9;

  return t1*(0.03109219370394*
          log(2.146079601698033*(0.00864486614818*
				 t7*t8*t10+0.0634682060977*t5*t6)/
	      (0.01855257090002*t7*t8/pow(t9,2.0)
	       +0.13620782246265*t5*t6*t10+1.0)+1.0)-0.062182*t3*t4);
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
  static num pbe_correlation(const num &R,
			     const num &S,
			     const num &Z,
			     const num &rhoa,
			     const num &rhob,
			     const num &zeta)
  {
    const num &sigma = Z;
    num d = 1.0/12*(sqrt(sigma)*pow(3,5.0/6))/
      (u(zeta)*pow(M_PI,-1.0/6)*pow(R,7.0/6));
    num eps = pw92_eps_pbe(R,zeta);
    return R*(eps + H(d,eps,zeta));
  }
}

using PBEc_internal::pbe_correlation;
using PBEx_internal::pbe_exchange;
using PBEx_internal::revpbe_exchange;


