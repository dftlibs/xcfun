#include "functional.h"
#include "slater.h"

template<class num>
static num ENERGY_FUNCTION(XC_SLATERX)(const densvars<num> &d)
{
    return slaterx(d); 
}      
NEW_LDA_FUNCTIONAL(XC_SLATERX);
SHORT_DESCRIPTION(XC_SLATERX) = "Slater LDA exchange";
LONG_DESCRIPTION(XC_SLATERX) = "LDA Exchange functional\n"
	     "P.A.M. Dirac, Proceedings of the Cambridge Philosophical "
	     "Society, 26 (1930) 376.\n"
	     "F. Bloch, Zeitschrift für Physik, 57 (1929) 545.\n\n"
	     "Implemented by Ulf Ekstrom\n"
	     "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n";
TEST_MODE(XC_SLATERX) = XC_PARTIAL_DERIVATIVES;
TEST_VARS(XC_SLATERX) = XC_A_B;
TEST_ORDER(XC_SLATERX) = 2;
TEST_THRESHOLD(XC_SLATERX) = 1e-11;
TEST_IN(XC_SLATERX) =     {0.39E+02, 0.38E+02};
TEST_OUT(XC_SLATERX) =    {-0.241948147838E+03, // energy
     -0.420747936684E+01, // gradient
     -0.417120618800E+01,
     -0.359613621097E-01, // hessian
     0,
     -0.365895279649E-01 };

