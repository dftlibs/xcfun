#include "functional.h"

  /* Short-range spin-dependent LDA exchange functional
     obtained by spin-scaling Ex[na,nb] = 1/2 (Ex[na]+Ex[nb])
     where Ex[n] is from Toulouse, Savin, Flad, IJQC 100, 1047 (2004)
     
     subroutine adapted from Molpro
     Created: 17-08-09, J. Toulouse, C++ by Ulf Ekstrom */

template<class num>
static num esrx_ldaerfspin(const num &na, parameter mu)
{

  const parameter ckf = 3.093667726280136;
  const num &rhoa = 2*na; // spin-scaling
  num akf = ckf*pow(rhoa,1.0/3.0);
  num a = mu/(2*akf);
  num a2 = a*a;
  num a3 = a2*a;
  // Test on the value of a. NB: These tests are not inifinitely
  // differentiable! 
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
static num ENERGY_FUNCTION(XC_LDAERFX)(const densvars<num> &d)
{
  double mu = d.get_param(XC_RANGESEP_MU);
  return 0.5*(esrx_ldaerfspin(d.a,mu) + esrx_ldaerfspin(d.b,mu)); 
}

NEW_LDA_FUNCTIONAL(XC_LDAERFX);
SHORT_DESCRIPTION(XC_LDAERFX) = "Short-range spin-dependent LDA exchange functional";
LONG_DESCRIPTION(XC_LDAERFX) =
	     "Short-range spin-dependent LDA exchange functional\n"
	     "obtained by spin-scaling Ex[na,nb] = 1/2 (Ex[na]+Ex[nb])\n"
	     "where Ex[n] is from Toulouse, Savin, Flad, IJQC 100, 1047 (2004)"
	     "Adapted from Gori-Giorgi and MOLPRO by Ulf Ekstrom\n"
	     "Test case from Gori-Giorgi (personal communication)\n"
	     "Range separation parameter is XC_RANGESEP_MU\n";
TEST_VARS(XC_LDAERFX) = XC_A_B;
TEST_MODE(XC_LDAERFX) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_LDAERFX) = 2;
TEST_THRESHOLD(XC_LDAERFX) = 1e-7;
TEST_IN(XC_LDAERFX) = {1.1, 1.0};
TEST_OUT(XC_LDAERFX) =    {
	-1.553573128702155,
	-1.067732841218878,
	-1.028091463003927,
	-0.3842706760115777,
	0,
	-0.4092115557248567	
    };

