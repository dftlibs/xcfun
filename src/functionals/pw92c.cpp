#include "functional.h"
#include "pw92eps.h"

template<class num>
static num ENERGY_FUNCTION(XC_PW92C)(const densvars<num> &d) 
{ 
  return pw92eps::pw92eps(d)*d.n;
}  

NEW_LDA_FUNCTIONAL(XC_PW92C);
SHORT_DESCRIPTION(XC_PW92C) = "PW92 LDA correlation";
LONG_DESCRIPTION(XC_PW92C) =	     "Electron-gas correlation energy\n"
	     "J.P.Perdew,Y. Wang; Phys. Rew. B; 40, 13244, (1992)\n"
	     "Implemented by Ulf Ekstrom. Some parameters have higher\n"
	     "accuracy than given in the paper.\n";
TEST_VARS(XC_PW92C) = XC_A_B;
TEST_MODE(XC_PW92C) = XC_PARTIAL_DERIVATIVES;
TEST_ORDER(XC_PW92C) = 2;
TEST_THRESHOLD(XC_PW92C) = 1e-11;
TEST_IN(XC_PW92C) = {0.39E+02, 0.38E+02};
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
TEST_OUT(XC_PW92C) =    {
      -8.4713855882783946e+00,
      -1.1861930857502517e-01,
      -1.2041769989725633e-01,
      +7.5202855619095870e-04,
      -1.0249091426230799e-03,	
      +7.9516089195232130e-04,
    };

