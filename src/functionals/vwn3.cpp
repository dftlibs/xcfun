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
#include "vwn.hpp"

template <typename num> static num vwn3c(const densvars<num> & d) {
  return d.n * vwn::vwn3_eps(d);
}

FUNCTIONAL(XC_VWN3C) = {
    "VWN3 LDA Correlation functional",
    "VWN3 LDA Correlation functional\n"
    "S.H. Vosko, L. Wilk, and M. Nusair: Accurate spin-dependent\n"
    "electron liquid correlation energies for local spin density\n"
    "calculations: a critical analysis, Can. J. Phys. 58 (1980) 1200-1211.\n"
    "Originally from Dalton, polished and converted by Ulf Ekstrom.\n",
    XC_DENSITY,
    ENERGY_FUNCTION(vwn3c)};
