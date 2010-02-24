
/*  Implementation of LYP functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
    autodiff version by Ulf Ekstrom 2009 */
template<>
struct functional<XC_LYP_CORRELATION>
{
  static const char *get_name(void) { return "lypc"; }
  static const char *get_reference(void) { return "blah blah"; }
  enum { supported_modes = XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
    { 
      const double A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
      using xc_constants::CF;
      
      num rho2 = d.n*d.n;
      num rhom13 = pow(d.n,-1.0/3.0);
      num denom = 1+D*rhom13;
      num omega = exp(-C*rhom13)/denom*pow(d.n,-11.0/3.0);
      num delta = rhom13*(C + D/denom);
      
      num t1 = pow(2.0,11.0/3.0)*CF*(pow(d.a,8.0/3.0) +
				     pow(d.b,8.0/3.0));
      t1 +=  ((47.0 - 7.0*delta)/18.0)*d.gnn;
      t1 += -(2.5 -delta/18.0)*(d.gaa+d.gbb);
      t1 +=  (11.0-delta)/9.0*(d.a*d.gaa + d.b*d.gbb)/d.n;
      num t5 = -2.0/3.0*rho2*d.gnn + 
	((2.0/3.0*rho2-d.a*d.a)*d.gbb +
	 (2.0/3.0*rho2-d.b*d.b)*d.gaa);
      return -A*(4*d.a*d.b/(denom*d.n)
		 +B*omega*(d.a*d.b*t1+t5));
    }
  static int test(void) 
  {
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_lyp.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      {
	-0.402158795173E+01,
	-0.762734644914E-01,
	-0.830226435821E-01,
	0.301052145436E-06,
	0.369624286402E-06,
	0.220298633297E-06,
	0.331769729999E-02,
	-0.248438749270E-02,
	-0.398359773843E-07,
	0.263970784129E-07,
	-0.335415277613E-08,
	0.384280348438E-02,
	0.275886078235E-07,
	-0.433118929134E-07,
	-0.685474898360E-08
      };
    return standard_abgga_test<XC_LYP_CORRELATION>(d,ref,1e-11);
  }
};

template<>
struct functional<XC_BECKE_EXCHANGE>
{
  static const char *get_name(void) { return "b88x"; }
  static const char *get_reference(void) { return "blah blah"; }
  enum { supported_modes = XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
    { 
      const double BECKE_THRESHOLD = 1e-14;
      const double BETA = 0.0042;  
      num ea = 0, eb = 0;
      
      if (d.b>BECKE_THRESHOLD)
	{
	  num xb = sqrt(d.gbb)*pow(d.b,-4.0/3.0);
	  num rb = pow(d.b,4.0/3.0);
	  num denomb = 1.0 +6.0*xb*BETA*asinh(xb);
	  eb = rb*xb*xb/denomb; 
	} 
      if (d.a>BECKE_THRESHOLD) 
	{
	  num xa = sqrt(d.gaa)*pow(d.a,-4.0/3.0);
	  num ra = pow(d.a,4.0/3.0);
	  num denoma = 1.0 +6.0*BETA*xa*asinh(xa);
	  ea = ra*xa*xa/denoma;
	}
      return -BETA*(ea+eb) + functional<XC_SLATER_EXCHANGE>::energy(d);
    }
  static int test(void) 
  {
    // Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_b88.html
    static const double d[5] = 
      {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
    static const double ref[21] =
      {	-0.277987329958E+03,
	-0.385951846654E+01,
	-0.381309494319E+01,
	-0.172434478018E-04,
	-0.173712338362E-04,
	0.000000000000E+00,
	-0.441426807406E-01,
	0.000000000000E+00,
	0.201415922856E-06,
	0.000000000000E+00,
	0.000000000000E+00,
	-0.447245742260E-01,
	0.000000000000E+00,
	0.195961359539E-06,
	0.000000000000E+00,
	0.700742719647E-11,
	0.000000000000E+00,
	0.000000000000E+00,
	0.718678968862E-11,
	0.000000000000E+00,
	0.000000000000E+00
      };
    return standard_abgga_test<XC_BECKE_EXCHANGE>(d,ref,1e-11);
  }
};

