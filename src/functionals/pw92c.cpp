#include "functional.h"
#include "pw92eps.h"

template<class num>
static num energy(const densvars<num> &d) 
{ 
  return pw92eps::pw92eps(d)*d.n;
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
  SET_GGA_ENERGY_FUNCTION(f,energy);
  const double d[] = 
    {0.39E+02, 0.38E+02};
  /* Daresbury numbers: 
  const double ref[] =
    {
      -0.847142825874E+01,
      -0.118619938062E+00,
      -0.120418335387E+00,
      0.752030427237E-03, 
      -0.102491320671E-02,
      0.795162900251E-03,
    };
  */
  /* Self computed numbers: */
  const double ref[] =
    {
      -8.4713855882783946e+00,
      -1.1861930857502517e-01,
      -1.2041769989725633e-01,
      +7.5202855619095870e-04,
      -1.0249091426230799e-03,	
      +7.9516089195232130e-04,
    };
  f.add_test(XC_VARS_AB,2,d,ref,1e-11);
}
