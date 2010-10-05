#ifndef TPSSX_EPS_H
#define TPSSX_EPS_H

#include "functional.h"
#include "constants.h"
#include "pbec_eps.h"

namespace tpssx_eps
{

  template<class num>
  static num epsx_unif(const densvars<num> &d)
  {
    return -3*cbrt(3*PI2*d.n)/(4*PI);
  }

  template<class num>
  static num x(const densvars<num> &d) 
  {    
    parameter kapa = 0.804;
    parameter mu = 0.21951;
    parameter b = 0.40;
    parameter e = 1.573;
    parameter c = 1.59096;
    num p = d.gnn/(4*pow(3*PI2,2.0/3.0)*pow(d.n,8.0/3.0));
    num z = d.gnn/(8.0*d.n*d.tau);
    num z2 = pow2(z);
//    tau_unif = 3*pow(3*PI2,2.0/3.0)*pow(d.n,5.0/3.0)/10.0;
    num alpha = 5*p*(1/z - 1)/3.0;
    num q_b = 9*(alpha-1)/(20*sqrt(1+ b*alpha*(alpha-1))) + 2*p/3.0;
    num x_a = p*(10.0/81.0 + c*z2/pow2(1+z2));
    x_a += 146.0*pow2(q_b)/2025.0;
    x_a -= 73.0*q_b*sqrt(0.5*0.6*0.6*z2 + 0.5*p*p);
    x_a += pow2(10.0*p/81.0)/kapa;
    x_a += sqrt(e)*0.6*0.6*z2*20.0/81.0 + e*mu*pow3(p);
    return x_a/(1+sqrt(e)*p);
  }

  template<class num>
  static num F_x(const densvars<num> &d)
  {
    parameter kapa = 0.804;
    num xpz = x(d);
    return 1 + kapa - kapa/(1+xpz/kapa);
  }

  template<class num>
  static num tpssx_eps(const densvars<num> &d)
  {
    num Fx = F_x(d);
    num epsxunif = epsx_unif(d);
    return epsxunif*Fx;
  }


}

#endif
