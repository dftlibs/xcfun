#include "functional.h"
#include "slater.h"

void setup_slaterx(functional &f)
{
  f.describe("slaterx",XC_LDA,
	     "Slater LDA exchange",
	     "LDA Exchange functional\n"
	     "P.A.M. Dirac, Proceedings of the Cambridge Philosophical "
	     "Society, 26 (1930) 376.\n"
	     "F. Bloch, Zeitschrift für Physik, 57 (1929) 545.\n\n"
	     "Implemented by Ulf Ekstrom\n"
	     "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n");

  SET_LDA_ENERGY_FUNCTION(f,slaterx);
  SET_GGA_ENERGY_FUNCTION(f,slaterx);

  static const double d[] = 
    {0.39E+02, 0.38E+02};
  static const double ref[] =
    {-0.241948147838E+03, // energy
     -0.420747936684E+01, // gradient
     -0.417120618800E+01,
     -0.359613621097E-01, // hessian
     0,
     -0.365895279649E-01 };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}
