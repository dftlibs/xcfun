#include "functional.h"
#include "pw9xx.h"

template<class num>
static num ENERGY_FUNCTION(XC_PW91X)(const densvars<num> &d)
{
  const parameter param_AB[6] = 
    { 0.19645, 7.7956, 0.2743, 0.15084, 100.0, 0.004};
  using pw91_like_x_internal::prefactor;
  using pw91_like_x_internal::pw91xk_enhancement;
  return prefactor(d.a)*pw91xk_enhancement(param_AB,d.a,d.gaa)
    + prefactor(d.b)*pw91xk_enhancement(param_AB,d.b,d.gbb);
}

NEW_GGA_FUNCTIONAL(XC_PW91X);
SHORT_DESCRIPTION(XC_PW91X) = "Perdew-Wang 1991 GGA Exchange Functional";
LONG_DESCRIPTION(XC_PW91X) =
	     "Perdew-Wang 1991 GGA Exchange Functional\n"
	     "J. P. Perdew, J. A. Chevary, S. H. Vosko, "
	     "K. A. Jackson, M. R. Pederson, and C. Fiolhais, "
	     "Phys. Rev. B 46, 6671 (1992)\n"
	     "Implemented by Andre Gomes.\n"
	     "Test from http://www.cse.scitech.ac.uk/ccg/dft/"
	     "data_pt_x_pw91.html\n"; 
TEST_VARS(XC_PW91X) = XC_A_B_GAA_GAB_GBB;
TEST_ORDER(XC_PW91X) = 2;
TEST_THRESHOLD(XC_PW91X) = 1e-11;
TEST_IN(XC_PW91X) =  {0.82E+02, 0.81E+02, 0.49E+07, 0.49E+07, 0.49E+07};
TEST_OUT(XC_PW91X) =    {
      -0.739934270280E+03,
      -0.500194130392E+01,
      -0.497593413511E+01,
      -0.661655297347E-05,
      0.000000000000E+00,
      -0.665149614704E-05,
      -0.259426653786E-01,
      0.000000000000E+00,
      0.352029178373E-07,
      0.000000000000E+00,
      0.000000000000E+00,
      -0.260706018375E-01,
      0.000000000000E+00,
      0.000000000000E+00,
      0.346740334540E-07,
      0.454242196579E-12,
      0.000000000000E+00,
      0.000000000000E+00,
      0.000000000000E+00,
      0.000000000000E+00,
      0.463780470889E-12,};

