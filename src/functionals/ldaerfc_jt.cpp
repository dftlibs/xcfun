#include "functional.hpp"
#include "vwn.hpp"

// This code was written by Radovan Bast based on the Ph.D. thesis of Julien Toulouse (page 159)
// and the Dalton implementation of the short-range spin-unpolarized LDA correlation functional
// by Julien Toulouse.

template<class num>
static num c1(const num &rs)
{
  const parameter u1 =  1.0270741452992294;
  const parameter u2 = -0.230160617208092;
  const parameter v1 =  0.6196884832404359;

  num rs2 = rs*rs;

  return (u1*rs + u2*rs2)/(1.0 + v1*rs);
}

template<class num>
static num c2(const densvars<num> &d)
{
  const parameter a   = 3.2581;
  const parameter f   = 3.39530545262710070631;
  const parameter bet = 163.44;
  const parameter gam = 4.7125;
  num   g0            = f*(pow(gam + d.r_s, 1.5) + bet)*exp(-a*sqrt(gam + d.r_s));
  num   n2            = d.n*d.n;
  num   denominator   = 0.5*M_PI*n2*(g0 - 0.5);
  num   result        = d.n*vwn::vwn5_eps(d)/denominator;

  return result;
}

template<class num>
static num ldaerfc_jt(const densvars<num> &d)
{
  double mu          = d.get_param(XC_RANGESEP_MU);
  num    denominator = 1.0 + c1(d.r_s)*mu + c2(d)*mu*mu;
  num    result      = d.n*vwn::vwn5_eps(d)/denominator;

  return result;
}

FUNCTIONAL(XC_LDAERFC_JT) = {
  "Short-range spin-unpolarized LDA correlation functional",
  "Short-range spin-unpolarized LDA correlation functional of\n"
  "Julien Toulouse et al.\n"
  "Written by Radovan Bast based on the Ph.D. thesis of Julien Toulouse (page 159)\n"
  "and the Dalton implementation by Julien Toulouse.\n"
  "Range separation parameter is XC_RANGESEP_MU\n",
  XC_DENSITY,
  ENERGY_FUNCTION(ldaerfc_jt)
};

// radovan:
// selftest yet to be written. i have compared SCF energy with Dalton
// and could obtain 11 matching digits


