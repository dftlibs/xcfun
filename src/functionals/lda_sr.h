#ifdef WITH_LDA_SR
/*
  This code was adapted by Ulf Ekstrom from a Fortran program provided
  by Gori-Giorgi.
 */

namespace ldaerf 
{

template<class num>
num esrc_ldaerfspin(const densvars<num> &d, double mu)
{
/*
C     Short-range spin-dependent LDA correlation functional from 
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
*/
  num eps = pw92_eps(d);
  return d.n*(eps - ecorrlr(d,mu,eps));
}

template<class num>
num Qrpa(const num &x)
{
  const double Acoul = 2.0*(log(2.0)-1.0)/(M_PI*M_PI);
  const double a2    = 5.84605;
  const double c2    = 3.91744;
  const double d2    = 3.44851;
  const double b2 = d2 - 3.0/(2*M_PI*Acoul)*pow(4.0/(9.0*M_PI),1.0/3.0);
  return Acoul*log( (1+x*(a2+x*(b2+c2*x)))/(1+x*(a2+d2*x)) );
}

template<class num>
num dpol(const num &rs)
{
  const double cf  = pow(9.0*M_PI/4.0,1.0/3.0);
  const double p2p = 0.04;
  const double p3p = 0.4319;
  num rs2 = rs*rs;
  return pow(2.0,5.0/3.0)/5.0*pow(cf,2)/rs2*
    (1.0 + (p3p - 0.454555)*rs)/(1.0 + p3p*rs + p2p*rs2);
}
/*
 on-top pair-distribution function
 Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
 x -> rs
*/
template<class num>
num g0f(const num &x)
{
  const double C0f = 0.0819306;
  const double D0f = 0.752411;
  const double E0f = -0.0127713;
  const double F0f = 0.00185898;
  return (1+ x*(D0f-0.7317 + x*(C0f + x*(E0f + F0f*x))))
	  *exp(-D0f*x)/2.0;
}

template<class num>
  num ecorrlr(const densvars<num> &d, double mu, const num &ec)
{
  const double alpha=pow(4.0/9.0/M_PI,1.0/3.0);
  double cf=1/alpha; //TODO: The normal CF?
  num phi=(pow(1.0+d.zeta,2.0/3.0)+pow(1.0-d.zeta,2.0/3.0))/2.0;
  // cc parameters from the fit
  const double adib = 0.784949;
  const double q1a = -0.388;
  const double q2a = 0.676;
  const double q3a = 0.547;
  const double t1a = -4.95;
  const double t2a = 1.0;
  const double t3a = 0.31;
  
  num b0 = adib*d.r_s;
  num rs2 = d.r_s*d.r_s;
  num rs3 = rs2*d.r_s;

  num d2anti=(q1a*d.r_s+q2a*rs2)*exp(-q3a*d.r_s)/rs2;
  num d3anti=(t1a*d.r_s+t2a*rs2)*exp(-t3a*d.r_s)/rs3;

  const num &z = d.zeta;
  num z2 = d.zeta*d.zeta;
  num coe2 = -3.0/8.0/rs3*(1.0-z2)*(g0f(d.r_s)-0.5);

  num coe3 = -(1.0-z2)*g0f(d.r_s)/(sqrt(2.0*M_PI)*rs3);

  num coe4 = -9.0/64.0/rs3*(pow((1.0+z)/2.0,2)*
			    dpol(d.r_s*pow(2/(1.0+z),1.0/3.0))+
			    pow((1.0-z)/2.0,2)*
			    dpol(d.r_s*pow(2.0/(1.0-z),1.0/3.0))+
			    (1-z*z)*d2anti-
			    pow(cf,2)/10.0*(pow(1.0+z,8.0/3.0)
					    +pow(1-z,8.0/3.0))/rs2);
  num coe5 = -9.0/40.0/(sqrt(2.0*M_PI)*rs3)*
    (pow((1.0+z)/2.0,2)*
     dpol(d.r_s*pow(2.0/(1.0+z),1.0/3.0))+
     pow((1.0-z)/2.0,2)*
     dpol(d.r_s*pow(2.0/(1.0-z),1.0/3.0))+(1.0-z2)*
     d3anti);

  num b06 = pow(b0,6);
  num b08 = pow(b0,8);
  num a1 = 4*b06*coe3 + b08*coe5;
  num a2 = 4*b06*coe2 + b08*coe4 + 6*pow(b0,4)*ec;
  num a3 = b08*coe3;
  num a4 = b06*(pow(b0,2)*coe2 + 4*ec);
    
  return (pow(phi,3)*Qrpa(mu*sqrt(d.r_s)/phi) + a1*pow(mu,3) + 
	  a2*pow(mu,4) + a3*pow(mu,5) +
          a4*pow(mu,6) + pow(b0*mu,8)*ec)/pow(1 + pow(b0*mu,2),4);
}


template<class num>
num esrx_ldaerfspin(const num &na, double mu)
{
/*
C     Short-range spin-dependent LDA exchange functional
C       obtained by spin-scaling Ex[na,nb] = 1/2 (Ex[na]+Ex[nb])
C       where Ex[n] is from Toulouse, Savin, Flad, IJQC 100, 1047 (2004)
C
C     subroutine adapted from Molpro
C     Created: 17-08-09, J. Toulouse, C++ by Ulf Ekstrom
*/
  const double ckf = 3.093667726280136;
  const num &rhoa = 2*na; // spin-scaling
  num akf = ckf*pow(rhoa,1.0/3.0);
  num a = mu/(2*akf);
  num a2 = a*a;
  num a3 = a2*a;
  // Test on the value of a. NB: These tests are not inifinitely differentiable
  if (a < 1e-9)
    // Limit for small a (expansion not so important as for large a)
    return -3.0/8.0*rhoa*pow(24.0*rhoa/M_PI,1.0/3.0);
  else if (a < 100)
    // Intermediate values of a
    return -(rhoa*pow(24.0*rhoa/M_PI,1.0/3.0))*
      (3.0/8.0 - a*(sqrt(M_PI)*erf(0.5/a) +
		    (2*a - 4*a3)*exp(-0.25/a2)-3.0*a + 4*a3));
  else if (a < 1e9)
    // Expansion for large a
    return - (rhoa*pow(24.0*rhoa/M_PI,1.0/3.0))/(96.0*a2);
  else
    // Limit for large a
    return 0;
}

template<class num>
num ldaerfspin_exchange(const densvars<num> &d, double mu)
{
  return 0.5*(esrx_ldaerfspin(d.a,mu) + esrx_ldaerfspin(d.b,mu)); 
}

}

template<>
struct functional<XC_SRLDA_ERF_CORRELATION>
{
  static const char *get_name(void) { return "srLDAerfc"; }
  static const char *get_reference(void) 
  { 
    return 
      "Short-range spin-dependent LDA correlation functional from\n"
      "Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)";
  }
  enum { supported_modes = XC_ALL_LDA | XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return ldaerf
      ::esrc_ldaerfspin(d, xc_settings[XC_SETTING_SRLDA_ERF_MU]);
  }
  static int test(void) 
  { 
    // Values provided by Paola Gori-Giori
    static const double test_mu = 0.4;
    static const double d[5] = {1.1, 1.0};
    static const double ref[21] = 
      {
	-0.1457945300494694,
	-0.07762517247521351,
	-0.08213242046922536,
	0,0,0,
	0.015795038615569225,
	-0.02744102459091926,
	0,0,0,
	0.019539653410626807};
    double mu_save = xc_settings[XC_SETTING_SRLDA_ERF_MU];
    xc_settings[XC_SETTING_SRLDA_ERF_MU] = test_mu;
    // Because of slightly different pw92 parameters in the reference
    // the results can only be compared to 7 digits.
    int res = standard_abgga_test<XC_SRLDA_ERF_CORRELATION>(d,ref,1e-7);
    xc_settings[XC_SETTING_SRLDA_ERF_MU] = mu_save;
    return res;
  }
};


template<>
struct functional<XC_SRLDA_ERF_EXCHANGE>
{
  static const char *get_name(void) { return "srLDAerfx"; }
  static const char *get_reference(void) 
  { 
    return 
      "Short-range spin-dependent LDA exchange functional\n"
      "obtained by spin-scaling Ex[na,nb] = 1/2 (Ex[na]+Ex[nb])\n"
      "where Ex[n] is from Toulouse, Savin, Flad, IJQC 100, 1047 (2004)"; 
  }
  enum { supported_modes = XC_ALL_LDA | XC_ALL_GGA };
  enum { max_order = XC_MAX_ORDER };
  template<class num>
  static num energy(const densvars<num> &d) 
  { 
    return ldaerf::
      ldaerfspin_exchange(d, xc_settings[XC_SETTING_SRLDA_ERF_MU]);
  }
  static int test(void) 
  {
    static const double test_mu = 0.4;
    static const double d[5] = {1.1, 1.0};
    static const double ref[21] = 
      {
	-1.553573128702155,
	-1.067732841218878,
	-1.028091463003927,
	0,0,0,
	-0.3842706760115777,
	0,
	0,0,0,
	-0.4092115557248567	
      };
    double mu_save = xc_settings[XC_SETTING_SRLDA_ERF_MU];
    xc_settings[XC_SETTING_SRLDA_ERF_MU] = test_mu;
    // Because of slightly different pw92 parameters in the reference
    // the results can only be compared to 7 digits.
    int res = standard_abgga_test<XC_SRLDA_ERF_EXCHANGE>(d,ref,1e-7);
    xc_settings[XC_SETTING_SRLDA_ERF_MU] = mu_save;
    return res;
  }
};
#endif
