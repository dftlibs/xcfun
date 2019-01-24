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

#include "functional.hpp"

template <typename num> static num vW_alpha(const num & na, const num & gaa) {
  return gaa / (8 * na);
}

template <typename num> static num vW(const densvars<num> & d) {
  return vW_alpha(d.a, d.gaa) + vW_alpha(d.b, d.gbb);
}

FUNCTIONAL(XC_VWK) = {"von Weizsaecker kinetic energy",
                      "von Weizsaecker kinetic energy\n"
                      "Implemented by Borgoo/Ekstrom.\n",
                      XC_DENSITY | XC_GRADIENT,
                      ENERGY_FUNCTION(vW)};
