/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2019 Ulf Ekström and contributors.
 *
 * This file is part of XCFun.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * XCFun library, see: <https://xcfun.readthedocs.io/>
 */

#include "pbex.hpp"
#include "functional.hpp"

// original PBE exchange functional

template <typename num> static num pbex_en(const densvars<num> & d) {
  return pbex::energy_pbe_ab(pbex::R_pbe, d.a, d.gaa) +
         pbex::energy_pbe_ab(pbex::R_pbe, d.b, d.gbb);
}

FUNCTIONAL(XC_PBEX) = {
    "PBE Exchange Functional",
    "PBE Exchange Functional\n"
    "J. P. Perdew, K. Burke, and M. Ernzerhof, "
    "Phys. Rev. Lett 77, 3865 (1996)\n"
    "Implemented by Ulf Ekstrom and Andre Gomes.\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(pbex_en)
#ifdef XCFUN_REF_PBEX_MU
        XC_A_B_GAA_GAB_GBB,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-11,
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06, 0.82E+06},
    {-0.276589791995E+03,

     -0.382556082420E+01, -0.378108116179E+01, -0.174145337536E-04,
     0.000000000000E+00,  -0.175120610339E-04,

     -0.429564214817E-01, 0.000000000000E+00,  0.185237729809E-06,
     0.000000000000E+00,  0.000000000000E+00,  -0.424802511645E-01,
     0.000000000000E+00,  0.000000000000E+00,  0.161839553501E-06,
     0.740514207206E-11,  0.000000000000E+00,  0.000000000000E+00,
     0.000000000000E+00,  0.000000000000E+00,  0.786563034093E-11}
#endif
};
