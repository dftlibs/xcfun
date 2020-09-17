/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2020 Ulf Ekström and contributors.
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

#include "SCAN_like_eps.hpp"
#include "constants.hpp"
#include "functional.hpp"

template <class num> static num r2SCANX(const densvars<num> & d) {
  num Fx_a;
  num Fx_b;
  Fx_a = SCAN_eps::get_SCAN_Fx(2.0 * d.a, 4.0 * d.gaa, 2.0 * d.taua, 2, 1, 1);
  num epsxunif_a = SCAN_eps::fx_unif(2 * d.a);
  Fx_b = SCAN_eps::get_SCAN_Fx(2.0 * d.b, 4.0 * d.gbb, 2.0 * d.taub, 2, 1, 1);
  num epsxunif_b = SCAN_eps::fx_unif(2 * d.b);

  return 0.5 * (Fx_a * epsxunif_a + Fx_b * epsxunif_b);
}

FUNCTIONAL(XC_R2SCANX) = {"r2SCAN exchange functional",
                          "r2SCAN exchange functional.\n"
                          "Restored-Regularised SCAN functional\n"
                          "J Furness, A Kaplan, J Ning, J Perdew, & J Sun; J. Chem. "
                          "Phys. Lett.; Accepted (DOI: 10.1021/acs.jpclett.0c02405)"
                          "Implemented by James Furness (@JFurness1)\n",
                          XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                          ENERGY_FUNCTION(r2SCANX) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                          XC_PARTIAL_DERIVATIVES,
                          1,
                          1e-11,
                          {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
                          {-1.62076202284e-01,
                           -8.60008389872e-01,
                           -5.57023311096e-01,
                           -3.79265043692e-02,
                           0.00000000000e+00,
                           -1.04606093722e-01,
                           5.89149529531e-02,
                           5.50319523444e-02}};
