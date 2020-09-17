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

template <class num> static num r4SCANC(const densvars<num> & d) {
  num eps_c = SCAN_eps::SCAN_C(d, 2, 1, 2);

  return eps_c;
}

FUNCTIONAL(XC_R4SCANC) = {"r4SCAN correlation functional",
                          "Equivalent to r2SCAN correlation functional.\n"
                          "Restored-Regularised fourt order SCAN functional\n"
                          "J Furness, in preparation"
                          "Implemented by James Furness (@JFurness1)\n",
                          XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                          ENERGY_FUNCTION(r4SCANC) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                          XC_PARTIAL_DERIVATIVES,
                          1,
                          1e-11,
                          {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
                          {-7.52536253477e-03,
                           -1.99923310008e-02,
                           -6.45082284064e-02,
                           1.10376734989e-02,
                           2.20753469977e-02,
                           1.10376734989e-02,
                           -1.76333423096e-02,
                           -1.76333423096e-02}};
