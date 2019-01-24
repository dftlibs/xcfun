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
#include "pbex.hpp"

template <typename num> static num rpbex(const densvars<num> & d) {
  using namespace pbex;
  using pw91_like_x_internal::prefactor;
  return prefactor(d.a) * enhancement_RPBE(d.a, d.gaa) +
         prefactor(d.b) * enhancement_RPBE(d.b, d.gbb);
}

FUNCTIONAL(XC_RPBEX) = {
    "RPBE Exchange Functional",
    "RPBE Exchange Functional\n"
    "Hammer, B. Hansen, L.B., Norskov, J.K.; PRB (59) p.7413, 1999\n"
    "Implemented by Ulf Ekstrom\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(rpbex)};
