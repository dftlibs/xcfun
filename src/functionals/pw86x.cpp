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

template <typename num> static num pw86x(const num & na, const num & gaa) {
  const parameter a = 1.0;
  const parameter b = 1.296;
  const parameter c = 14.0;
  const parameter d = 0.20;
  num rho = 2 * na, grad2 = 4 * gaa;
  const num Ax = -pow(3.0 / M_PI, 1.0 / 3.0) * 3.0 / 4.0;
  const num kf = pow(3.0 * pow(M_PI, 2) * rho, 1.0 / 3.0);
  num s2 = grad2 / pow(2.0 * kf * rho, 2);
  num F = pow(a + s2 * (b + s2 * (c + d * s2)), 1.0 / 15.0);
  return Ax * pow(rho, 4.0 / 3.0) * F;
}

template <typename num> static num pw86xtot(const densvars<num> & d) {
  return 0.5 * (pw86x(d.a, d.gaa) + pw86x(d.b, d.gbb));
}

FUNCTIONAL(XC_PW86X) = {
    "PW86 exchange",
    "Perdew-Wang 86 GGA exchange including Slater part\n"
    "Phys. Rev. B 33. 8800 (1986)\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(pw86xtot) XC_A_B_GAA_GAB_GBB,
};
