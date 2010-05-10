#include "functional.h"
#include "pw9xx.h"

template<class num>
static num energy(const densvars<num> &d)
{
  static const parameter param_AB[6] = 
    { 0.19645, 7.7956, 0.2743, 0.1508, 100.0, 0.004};
  using pw91_like_x_internal::prefactor;
  using pw91_like_x_internal::pw91xk_enhancement;
  return prefactor(d.a)*pw91xk_enhancement(param_AB,d.a,d.gaa)
    + prefactor(d.b)*pw91xk_enhancement(param_AB,d.b,d.gbb);
}

void setup_pw91x(functional &f)
{
  f.describe("pw91x",XC_GGA,
	     "Perdew-Wang 1991 GGA Exchange Functional",
	     "Perdew-Wang 1991 GGA Exchange Functional\n"
	     "J. P. Perdew, J. A. Chevary, S. H. Vosko, "
	     "K. A. Jackson, M. R. Pederson, and C. Fiolhais, "
	     "Phys. Rev. B 46, 6671 (1992)\n"
	     "Implemented by Andre Gomes.\n"
	     "Test from http://www.cse.scitech.ac.uk/ccg/dft/"
	     "data_pt_x_pw91.html\n"); 
  SET_GGA_ENERGY_FUNCTION(f,energy);
  static const double d[] = 
    {0.82E+02, 0.81E+02, 0.49E+07, 0.49E+07, 0.49E+07};
  static const double out[] = 
    {
      -0.739934270280E+03,
      -0.500194130392E+01,
      -0.497593413511E+01,
      -0.661655297347E-05,
      -0.665149614704E-05,
      0.000000000000E+00,
      -0.259426653786E-01,
      0.000000000000E+00,
      0.352029178373E-07,
      0.000000000000E+00,
      0.000000000000E+00,
      -0.260706018375E-01,
      0.000000000000E+00,
      0.346740334540E-07,
      0.000000000000E+00,
      0.454242196579E-12,
      0.000000000000E+00,
      0.000000000000E+00,
      0.463780470889E-12,
      0.000000000000E+00,
      0.000000000000E+00};
  f.add_test(XC_VARS_AB,2,d,out,1e-11);
}
