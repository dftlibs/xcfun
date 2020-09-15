#include "SCAN_like_eps.hpp"
#include "constants.hpp"
#include "functional.hpp"

template <class num> static num rppSCANC(const densvars<num> & d) {
  num eps_c = SCAN_eps::SCAN_C(d, 2, 1, 0);

  return eps_c;
}

FUNCTIONAL(XC_RPPSCANC) = {"r++SCAN correlation functional",
                           "r++SCAN correlation functional.\n"
                           "Restored-Regularised SCAN functional\n"
                           "J Furness, in preparation"
                           "Implemented by James Furness (@JFurness1)\n",
                           XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                           ENERGY_FUNCTION(rppSCANC) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                           XC_PARTIAL_DERIVATIVES,
                           1,
                           1e-11,
                           {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
                           {-7.52313999333e-03,
                            -1.98575241135e-02,
                            -6.43444019823e-02,
                            1.09821296406e-02,
                            2.19642592812e-02,
                            1.09821296406e-02,
                            -1.76156638566e-02,
                            -1.76156638566e-02}};
