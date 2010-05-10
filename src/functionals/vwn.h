#ifndef VWN_H
#define VWN_H

#include "functional.h"

namespace vwn
{

  static parameter vwn_a(const parameter p[])
  {
    return p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1;
  }
  
  static parameter vwn_b(const parameter p[])
  {
    return 2*(p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1) + 2;
  }
  
  static parameter vwn_c(const parameter p[])
  {
    return 2*p[2]*(1/sqrt(4*p[3] - p[2]*p[2]) - p[0]/
		   ((p[0]*p[0] + p[0]*p[2] + p[3])*
		    sqrt(4*p[3] - p[2]*p[2])/(p[2] + 2*p[0])));
  }
  
  template<class num>
  static num vwn_x(const num& s,const parameter p[])
  {
    return s*s + p[2]*s + p[3];
  }
  
  template<class num>
  static num vwn_y(const num& s,const parameter p[])
  {
    return s - p[0];
  }
  
  template<class num>
  static num vwn_z(const num& s,const parameter p[])
  {
    return sqrt(4*p[3] - p[2]*p[2])/(2*s + p[2]);
  }
  
  template<class num>
  static num vwn_f(const num& s,const parameter p[])
  {
    return 0.5*p[1]*(2*log(s) + vwn_a(p)*log(vwn_x(s, p)) - 
		     vwn_b(p)*log(vwn_y(s, p)) + 
		     vwn_c(p)*atan(vwn_z(s, p)));
  }

  template<class num>
  static num vwn5c(const densvars<num> &d)
  {
    //ulfek: second elements are multiplied by 2 wrt molpro manual
    static const parameter  para[] = {-0.10498,   0.0621814, 3.72744, 12.9352};
    static const parameter ferro[] = {-0.325,     0.0310907, 7.06042, 18.0578};
    static const parameter inter[] = {-0.0047584,-pow(3*M_PI*M_PI,-1.0), 1.13107, 13.0045};
  // ulfek: in Dalton and Dirac inter[1] has too few decimals, should
  // be -1/(3pi^2)
    num s = sqrt(d.r_s);
    // Constant is (2^1/3-1)^-1/2
    num g = 1.92366105093154*(pow(1 + d.zeta,4.0/3) + pow(1 - d.zeta,4.0/3) - 2);
    num zeta4 = pow(d.zeta,4);
    num dd = g*((vwn_f(s, ferro) - vwn_f(s, para))*zeta4 +
		vwn_f(s, inter)*(1 - zeta4)*(9.0/4.0*(pow(2,1.0/3.0)-1)));
    // Dalton (for some reason??) uses a value:    0.584822305543806
    // The real value should be 9/4 (2^(1/3) -1) = 0.5848223622634646
    return d.n*(vwn_f(s, para) + dd);
  }
}

#endif
