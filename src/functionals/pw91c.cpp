#include "functional.hpp"
#include "constants.hpp"


template<class num>
static num Gc(const num &r, 
	      const parameter A, 
	      const parameter a1,
	      const parameter b1,
	      const parameter b2,
	      const parameter b3,
	      const parameter b4,
	      const parameter p)
{
  num sqrtr = sqrt(r);
  return -2*A*(1+a1*r)*log(1 + 0.5/(A*(sqrtr*(b1 + sqrtr*(b2 + sqrtr*b3)) + b4*pow(r,p+1))));
}

template<class num>
static num pw91c(const densvars<num> &d) 
{
  const parameter pa = 1.0;
  const parameter Aa = 0.016887;
  const parameter a1a = 0.11125;
  const parameter b1a = 10.357;
  const parameter b2a = 3.6231;
  const parameter b3a = 0.88026;
  const parameter b4a = 0.49671;
  const parameter pe = 1;
  const parameter c0p = 0.031091;
  const parameter a1p = 0.21370;
  const parameter b1p = 7.5957;
  const parameter b2p = 3.5876;
  const parameter b3p = 1.6382;
  const parameter b4p = 0.49294;
  const parameter c0f = 0.015545;
  const parameter a1f = 0.20548;
  const parameter b1f = 14.1189;
  const parameter b2f = 6.1977;
  const parameter b3f = 3.3662;
  const parameter b4f = 0.62517;
  const parameter d2fz0 = 1.709921; 
  num fz  = (ufunc(d.zeta,4.0/3.0)-2)/(2*pow(2,1.0/3.0)-2);
  num Ac  = Gc(d.r_s,  Aa, a1a, b1a, b2a, b3a, b4a, pa);
  num EcP = Gc(d.r_s, c0p, a1p, b1p, b2p, b3p, b4p, pe);
  num EcF = Gc(d.r_s, c0f, a1f, b1f, b2f, b3f, b4f, pe);
  num Ec = EcP - Ac*fz*(1-pow(d.zeta,4))/d2fz0 + (EcF - EcP)*fz*pow(d.zeta,4);
  num kF = cbrt(3*M_PI*M_PI*d.n);
  num ks = sqrt(4)*sqrt(kF/M_PI);
  num gs = 0.5*ufunc(d.zeta,2.0/3.0);
  num T2 = 0.25*d.gnn/pow(gs*ks*d.n,2);
  const parameter alpha = 0.09;
  const parameter Cc0 = 0.004235;
  const parameter Cx = -0.001667;
  const parameter nu = 16*cbrt(3*M_PI*M_PI)/M_PI;
  const parameter beta = nu*Cc0;
  num A = 2*alpha/beta/expm1(-2*(alpha*Ec)/(pow(gs,3)*beta*beta));
  num Cc = 1.0/1000*((2.568 + d.n*(23.266 + 0.007389*d.n))/
		     (1 + d.n*(8.723 + d.n*(0.472 + d.n*0.073890)))) - Cx;
  num H0 = 0.5*pow(gs,3)*beta*beta*log((1 + 2*(alpha*(T2+A*T2*T2))/(beta*(1 + A*T2*(1 + A*T2)))))/alpha;
  num H1 = nu*(Cc - Cc0 - 3.0/7.0*Cx)*pow(gs,3)*T2*exp(-100*pow(gs,4)*pow(ks,2)*T2/pow(kF,2));  
  return d.n*(Ec + H0 + H1);
}

FUNCTIONAL(XC_PW91C) = {
  "PW91 Correlation",
  "J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, M.R. Pederson, D.J. Singh and"
  "C. Fiolhais, 'Atoms, molecules, solids and surfaces: Applications of the generalized"
  "gradient approximation for exchange and correlation', Phys. Rev. B, 46(11):6671"
  "6687, 1992"
  "Implemented by Ulf Ekstrom. Test from ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_c_pw91.html\n",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(pw91c)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2, // test order
  1e-11, // test threshold
  {0.78E-01,0.31E-01, 0.41E-02, 0.38E-02, 0.36E-02},
  {
    -0.450106022368E-02,
    -0.568426120793E-01,
    -0.876498238182E-01,
    0.531751834129E-01,
    0.106350366826E+00,
    0.531751834129E-01,
    0.192611490688E+00,
    -0.331153469676E+00,
    -0.308498888708E+00,
    -0.616997777417E+00,
    -0.308498888708E+00,
    0.645433404561E+00,
    -0.129198778272E+00,
    -0.258397556545E+00,
    -0.129198778272E+00,
    -0.119074508391E+01,
    -0.238149016782E+01,
    -0.119074508391E+01,
    -0.476298033563E+01,
    -0.238149016782E+01,
    -0.119074508391E+01,
  },
};
