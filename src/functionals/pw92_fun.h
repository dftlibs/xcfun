#ifndef PW92_H
#define PW92_H

namespace pw92_internal
{
  template<class num>
  static num e(const num &r, const double t[7])
  {
    return -2*t[0]*(1+t[1]*r)*
      log(1+0.5/(t[0]*(t[2]*sqrt(r)+t[3]*r+t[4]*
		       pow(r,3.0/2)+t[5]*pow(r,t[6]+1))));
  }

  template<class num>
  static num omega(const num &z)
  {
    return (pow(1+z,4.0/3)+pow(1-z,4.0/3)-2)/(2*pow(2,1.0/3)-2);
  }

  template<class num>
  static num pw92_eps(const densvars<num> &d)
  {
    // The T..P parameters, but transposed for easier use with the e function
    const double TUVWXYP[3][7] = 
      {{0.031091, 0.21370, 7.5957 , 3.5876, 1.6382 ,0.49294,1},
       {0.015545 , 0.20548, 14.1189, 6.1977, 3.3662 ,0.62517,1},
       {0.016887, 0.11125, 10.357 , 3.6231, 0.88026,0.49671,1}};
    const double c = 1.709921;
    num zeta4 = pow(d.zeta,4);
    num omegaval = omega(d.zeta);
    return e(d.r_s,TUVWXYP[0])
      -  e(d.r_s,TUVWXYP[2])*omegaval*(1-zeta4)/c 
      + (e(d.r_s,TUVWXYP[1]) - e(d.r_s,TUVWXYP[0]))*omegaval*zeta4;
  }
  // This is with the Daresbury parameters for PBEc, not exactly
  // the same as for Daresbury pw92. 
  template<class num>
  static num pw92_eps_pbe(const densvars<num> &d)
  {
    // The T..P parameters, but transposed for easier use with the e function
    const double TUVWXYP[3][7] = 
      {{0.0310907, 0.21370, 7.5957 , 3.5876, 1.6382 ,0.49294,1},
       {0.01554535 , 0.20548, 14.1189, 6.1977, 3.3662 ,0.62517,1},
       {0.0168869, 0.11125, 10.357 , 3.6231, 0.88026,0.49671,1}};
    //    const double c = 1.709921;
    const double c = 8.0/(9.0*(2*pow(2,1.0/3.0)-2));
    num zeta4 = pow(d.zeta,4);
    num omegaval = omega(d.zeta);
    return e(d.r_s,TUVWXYP[0])
      -  e(d.r_s,TUVWXYP[2])*omegaval*(1-zeta4)/c 
      + (e(d.r_s,TUVWXYP[1]) - e(d.r_s,TUVWXYP[0]))*omegaval*zeta4;
  }
}


template<class num>
static num pw92_correlation(const densvars<num> &d)
{
  return d.n*pw92_eps(d);
}

using pw92_internal::pw92_eps;
using pw92_internal::pw92_eps_pbe;

template<>
struct functional<XC_PW92_CORRELATION>
{
  static const char *get_name(void) { return "pw92c"; }
  static const char *get_reference(void) 
  {
    return "PW92 correlation functional\n"
      "Electron-gas correlation energy\n"
      "J.P.Perdew,Y. Wang; Phys. Rew. B; 40, 13244, (1992)";
; }
  enum { supported_modes = XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return pw92_internal::pw92_eps(d)*d.n;
  }

  static int test(void) 
  {
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_b88.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      {
	-0.847142825874E+01,
	-0.118619938062E+00,
	-0.120418335387E+00,
	0.000000000000E+00,
	0.000000000000E+00,
	0.000000000000E+00,
	0.752030427237E-03, // Suspect
	-0.102491320671E-02,// Suspect
	0,0,0,
	0.795162900251E-03,// Suspect
      };
    return standard_abgga_test<XC_PW92_CORRELATION>(d,ref,1e-11);
  }
};

#endif
