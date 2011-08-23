#include "functional.h"
#include "pz81c.h"

template<class num>
static num pz81c(const densvars<num> &d) 
{ 
  return pz81eps::pz81eps(d)*d.n;
}

FUNCTIONAL(XC_PZ81C) = {
  "PZ81 LDA correlation",  
  "Implemented by Ulf Ekstrom. Test from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_pz81.html\n",
  XC_DENSITY,
  ENERGY_FUNCTION(pz81c)
  XC_A_B,
  XC_PARTIAL_DERIVATIVES,
  2,
  1e-11,
  {0.39E+02, 0.38E+02},
  { 
    -0.847966726388E+01,
    -0.118689256817E+00,
    -0.121016447321E+00,
    0.100782705799E-02,
    -0.129204915717E-02,
    0.106274642401E-02,
  }
};

