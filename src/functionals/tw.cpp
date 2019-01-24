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

#include "constants.hpp"
#include "functional.hpp"

//  von Weizsacker kinetic energy functional

template <typename num> static num tw(const densvars<num> & d) {
  using xc_constants::CF;

  return 1. / 8. * pow(d.gaa + d.gbb, 2.0) / d.n;
}

FUNCTIONAL(XC_TW) = {
    "von Weizsacker Kinetic Energy Functional",
    "von Weizsacker Kinetic Energy Functional\n"
    "Implemented by AB and SR.\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(tw) XC_A_B_GAA_GAB_GBB,
};
