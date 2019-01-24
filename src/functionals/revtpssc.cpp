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
#include "revtpssc_eps.hpp"

template <typename num> static num revtpssc(const densvars<num> & d) {
  num eps = revtpssc_eps::revtpssc_eps(d);
  return d.n * eps;
}

FUNCTIONAL(XC_REVTPSSC) = {
    "Revised TPSS correlation functional",
    "Revised TPSS correlation functional.\n"
    "J.P. Perdew, A. Ruzsinszky, G. I. Csonka, L. A. Constantin, J. Sun,\n"
    "Workhorse Semilocal Density Functional\n"
    "for Condensed Matter Physics and Quantum Chemistry\n"
    "Phys. Rev. Lett. 103 (2009) 026403\n"
    "Implemented by Andrea Debnarova\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
    ENERGY_FUNCTION(revtpssc)};
