#include "functional.h"
#include "pw92eps.h"

template<class num>
static num energy(const densvars<num> &d) 
{ 
  return pw92eps::eps(d)*d.n;
}  

void setup_pw92c(functional &f)
{
  f.describe(XC_PW92C, XC_LDA,
	     "PW92 LDA correlation",
	     "Electron-gas correlation energy\n"
	     "J.P.Perdew,Y. Wang; Phys. Rew. B; 40, 13244, (1992)\n"
	     "Implemented by Ulf Ekstrom. Some parameters have higher\n"
	     "accuracy than given in the paper.\n");
  SET_LDA_ENERGY_FUNCTION(f,energy);
  static const double d[] = 
    {0.39E+02, 0.38E+02};
  static const double ref[] =
    {
      -0.847142825874E+01,
      -0.118619938062E+00,
      -0.120418335387E+00,
      0.752030427237E-03, 
      -0.102491320671E-02,
      0.795162900251E-03,
    };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}
