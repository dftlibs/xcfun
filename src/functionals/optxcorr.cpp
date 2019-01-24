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

template <typename num> static num optxcorr(const densvars<num> & d) {

  /*AMT Note this implements just the correction part of optx without the weighting
     parameters
        hence the sign change on the return value below relative to optx.cpp.
        This can be combined with specify slaterx, ktx and lypc to specify e.g. the
     KT3 functional
        as: 1.092*slaterx - 0.004*ktx - 0.925452*optxcorr + 0.864409*lypc
  */

  const parameter gamma = 0.006;
  num g_xa2 = gamma * d.gaa * pow(d.a, -8.0 / 3.0);
  num g_xb2 = gamma * d.gbb * pow(d.b, -8.0 / 3.0);
  return (d.a_43 * (pow(g_xa2, 2) * pow(1 + g_xa2, -2))) +
         (d.b_43 * (pow(g_xb2, 2) * pow(1 + g_xb2, -2)));
}

FUNCTIONAL(XC_OPTXCORR) = {
    "OPTX Handy & Cohen exchange -- correction part only",
    "OPTX Handy & Cohen exchange GGA exchange functional -- correction part only\n"
    "Implemented by AMT\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(optxcorr)};
