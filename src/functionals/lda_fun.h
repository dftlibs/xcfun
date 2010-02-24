
namespace VWN5_internal
{

static double vwn_a(const double p[])
{
  return p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1;
}

static double vwn_b(const double p[])
{
  return 2*(p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1) + 2;
}

static double vwn_c(const double p[])
{
  return 2*p[2]*(1/sqrt(4*p[3] - p[2]*p[2]) - p[0]/
		 ((p[0]*p[0] + p[0]*p[2] + p[3])*
		  sqrt(4*p[3] - p[2]*p[2])/(p[2] + 2*p[0])));
}

template<class num>
static num vwn_x(const num& s,const double p[])
{
  return s*s + p[2]*s + p[3];
}

template<class num>
static num vwn_y(const num& s,const double p[])
{
  return s - p[0];
}

template<class num>
static num vwn_z(const num& s,const double p[])
{
  return sqrt(4*p[3] - p[2]*p[2])/(2*s + p[2]);
}

template<class num>
static num vwn_f(const num& s,const double p[])
{
  return 0.5*p[1]*(2*log(s) + vwn_a(p)*log(vwn_x(s, p)) - 
		   vwn_b(p)*log(vwn_y(s, p)) + 
		   vwn_c(p)*atan(vwn_z(s, p)));
}

template<class num>
static num vwn5_correlation(const densvars<num> &d)
{
  //ulfek: second elements are multiplied by 2 wrt molpro manual
  static const double  para[] = {-0.10498,   0.0621814, 3.72744, 12.9352};
  static const double ferro[] = {-0.325,     0.0310907, 7.06042, 18.0578};
  static const double inter[] = {-0.0047584,-pow(3*M_PI*M_PI,-1.0), 1.13107, 13.0045};
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

template<>
struct functional<XC_VWN5_CORRELATION>
{
  static const char *get_name(void) { return "vwn5c"; }
  static const char *get_reference(void) { return "blah blah"; }
  enum { supported_modes = XC_ALL_LDA | XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  {
    return VWN5_internal::vwn5_correlation(d);
  }
  static int test(void) 
  {
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_vwn5.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      { -0.851077910672E+01,
	-0.119099058995E+00,
	-0.120906044904E+00,
	0.000000000000E+00,
	0.000000000000E+00,
	0.000000000000E+00,
	0.756836181702E-03,
	-0.102861281830E-02,
	0,0,0,
	0.800136175083E-03};
    return standard_abgga_test<XC_VWN5_CORRELATION>(d,ref,1e-11);
  }
};

template<>
struct functional<XC_SLATER_EXCHANGE>
{
  static const char *get_name(void) { return "slaterx"; }
  static const char *get_reference(void) { return "BLAH BLAH"; }
  enum { supported_modes = XC_ALL_LDA | XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return -3.0/4.0*pow(6/M_PI, 1.0/3.0)*
      (pow(d.a,4.0/3.0) + pow(d.b,4.0/3.0));
  }
  static int test(void) 
  { 
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      {-0.241948147838E+03, // energy
       -0.420747936684E+01, // gradient
       -0.417120618800E+01,
       0.000000000000E+00,
       0.000000000000E+00,
       0.000000000000E+00,
       -0.359613621097E-01, // hessian
       0.000000000000E+00,
       0,0,0,
      -0.365895279649E-01 }; // Rest are zero because this is LDA
    return standard_abgga_test<XC_SLATER_EXCHANGE>(d,ref,1e-11);
  }
};
