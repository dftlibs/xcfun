#include "functional.h"
#include "constants.h"

//  Thomas-Fermi kinetic energy functional

template<class num>
static num energy(const densvars<num> &d)
{
    using xc_constants::CF;

    return CF*pow(d.n, 5.0/3.0);
}

void setup_tfk(functional &f)
{
  f.describe("tfk",XC_LDA,
	     "Thomas-Fermi Kinetic Energy Functional",
	     "Thomas-Fermi Kinetic Energy Functional\n"
	     "\n"
	     "Implemented by Andre Gomes.\n"); 
  SET_LDA_ENERGY_FUNCTION(f,energy);
  static const double d[] =
    {1., .8};
  static const double ref[] =
    { 7.64755771625168, 7.08107195949229, 7.08107195949229, 2.62261924425641, 2.62261924425641, 2.62261924425641 };
  f.add_test(XC_VARS_AB,1,d,ref,1e-5);

}
