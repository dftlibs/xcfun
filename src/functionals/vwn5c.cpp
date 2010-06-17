#include "functional.h"
#include "vwn.h"

void setup_vwn5(functional &f)
{
  f.describe(XC_VWN5C, XC_LDA,
	     "VWN5 LDA Correlation functional",
	     "VWN5 LDA Correlation functional\n"
	     "S.H. Vosko, L. Wilk, and M. Nusair: Accurate spin-dependent\n"
	     "electron liquid correlation energies for local spin density\n"
	     "calculations: a critical analysis, Can. J. Phys. 58 (1980) 1200-1211.\n"
	     "Originally from Dalton, polished and converted by Ulf Ekstrom.\n"
	     "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_vwn5.html\n");

  SET_LDA_ENERGY_FUNCTION(f,vwn::vwn5c);
  SET_GGA_ENERGY_FUNCTION(f,vwn::vwn5c);

  static const double d[] = {0.39E+02, 0.38E+02};
  static const double ref[] =
    {
      -0.851077910672E+01,
      -0.119099058995E+00,
      -0.120906044904E+00,
      0.756836181702E-03,
      -0.102861281830E-02,
      0.800136175083E-03
    };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}
