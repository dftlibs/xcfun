#include "functional.hpp"
#include "constants.hpp"
#include "pz81c.hpp"

template<class num>
static num Cg(const num &r) 
{ 
  parameter Cx = 0.001667;
  parameter Bg = 0.000007389;
  return Cx + (0.002568 + r*(0.023266 + Bg*r))/(1 + r*(8.723 + r*(0.472 + 10000*Bg*r)));
}

template<class num>
static num Pg(const densvars<num> &d) 
{ 
  parameter Fg = 0.11;
  parameter Cinf = 0.004235;
  parameter fudge = 1e-12; // Avoid instability at d.gnn = 0
  return 1.745*Fg*Cinf*sqrt(fudge + d.gnn)/(Cg(d.r_s)*pow(d.n,7.0/6.0));
}

template<class num>
static num dz(const densvars<num> &d) 
{ 
  return cbrt(2.0)*sqrt(pow(d.a,5.0/3.0) + pow(d.b,5.0/3.0))*pow(d.n,-5.0/6.0);
}

template<class num>
static num p86c(const densvars<num> &d) 
{ 
  return d.n*pz81eps::pz81eps(d) + exp(-Pg(d))*Cg(d.r_s)*d.gnn/(pow(d.n,4.0/3.0)*dz(d));
}

FUNCTIONAL(XC_P86C) = {
  "P86C GGA correlation",  
  "J.P. Density-functional approximation for the correlation energy\n"
  "of the inhomogeneous electron , Phys. Rev. B, 33(12):8822gasPerdew,\n" 
  "Implemented by Ulf Ekstrom.\n"
  "Reference data from ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_c_p86.html",
  XC_DENSITY | XC_GRADIENT,
  ENERGY_FUNCTION(p86c)
  XC_A_B_GAA_GAB_GBB,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06},
  { 
    -0.356963343227E+01,
    -0.433530899660E-01,
    -0.447602011737E-01,
    -0.122488989075E-06,
    -0.244977978150E-06,
    -0.122488989075E-06,
    -0.169450157064E-02,
    -0.308422706942E-02,
    0.227151852191E-07,
    0.454303704383E-07,
    0.227151852191E-07,
    -0.165896811718E-02,
    0.226922323054E-07,
    0.453844646108E-07,
    0.226922323054E-07,
    -0.207808284633E-12,
    -0.415616569266E-12,
    -0.207808284633E-12,
    -0.831233138532E-12,
    -0.415616569266E-12,
    -0.207808284633E-12,
  }
};
