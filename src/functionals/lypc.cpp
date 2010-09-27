#include "functional.h"
#include "constants.h"

template<class num>
static num new_energy(const densvars<num> &d) 
{
  const parameter A = 0.04918;
  const parameter B = 0.132;
  const parameter C = 0.2533;
  const parameter Dd = 0.349;
  using xc_constants::CF;
  num icbrtn = pow(d.n,-1.0/3.0);
  num P = 1/(1+Dd*icbrtn);
  num omega = exp(-C*icbrtn)*P*pow(d.n,-11.0/3.0);
  num delta = icbrtn*(C+Dd*P);
  num n2 = d.n*d.n;
  return 
    -A*(4*d.a*d.b*P/d.n +
	B*omega*(d.a*d.b*(pow(2,11.0/3.0)*CF*(pow(d.a,8.0/3.0)+pow(d.b,8.0/3.0))
			  +(47.0 - 7.0*delta)*d.gnn/18.0
			  -(2.5 - delta/18.0)*(d.gaa + d.gbb)
			  -(delta-11.0)/9.0*(d.a*d.gaa + d.b*d.gbb)/d.n)
		 - 2.0/3.0*n2*d.gnn 
		 + (2.0/3.0*n2 - d.a*d.a)*d.gbb
		 + (2.0/3.0*n2 - d.b*d.b)*d.gaa));
}

void setup_lypc(functional &f)
{
  f.describe(XC_LYPC,XC_GGA,
	     "LYP correlation",
	     "LYP correlation\n"
	     "C. Lee, W. Yang, and R.G. Parr, Development of the \n"
	     "Colle-Salvetti correlation-energy formula into a functional\n"
	     "of the electron density, Phys. Rev. B37 (1988) 785-789\n"
	     "Implemented by Ulf Ekstrom\n"
	     "Test: http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_lyp.html\n");
  SET_GGA_ENERGY_FUNCTION(f,new_energy);
  const double d[5] = 
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06,0.82E+06};
  const double ref[21] =
    {
      -0.402158795173E+01,

      -0.762734644914E-01,
      -0.830226435821E-01,
      0.301052145436E-06,
      0.220298633297E-06,
      0.369624286402E-06,
      0.331769729999E-02,

      -0.248438749270E-02,
      -0.398359773843E-07,
      -0.335415277613E-08,
      0.263970784129E-07,
      0.384280348438E-02,
      0.275886078235E-07,
      -0.685474898360E-08,
      -0.433118929134E-07,
    };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}



