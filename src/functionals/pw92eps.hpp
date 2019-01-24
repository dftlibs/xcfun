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

#pragma once

#include "constants.hpp"

namespace pw92eps {
template <typename num> static num eopt(const num & sqrtr, const parameter t[6]) {
  return -2 * t[0] * (1 + t[1] * sqrtr * sqrtr) *
         log(1 + 0.5 / (t[0] *
                        (sqrtr *
                         (t[2] + sqrtr * (t[3] + sqrtr * (t[4] + t[5] * sqrtr))))));
}

template <typename num> static num omega(const densvars<num> & d) {
  return (cbrt(2.0) / (cbrt(2.0) - 1)) * (d.a_43 + d.b_43) * d.n_m13 / d.n -
         1 / (cbrt(2.0) - 1);
}

template <typename num> static num omega(const num & z) {
// return (pow(1+z,4.0/3.0)+pow(1-z,4.0/3.0)-2)/(2*pow(2,1.0/3.0)-2);
#ifndef XCFUN_REF_PW92C // has effect of about 10^-11 in the testcase energy
  return (ufunc(z, 4.0 / 3.0) - 2) / (2 * pow(2, 1.0 / 3.0) - 2);
#else
  return (ufunc(z, 4.0 / 3.0) - 2) / 0.5198421;
#endif
}

#define PW92C_PARAMS                                                                \
  {{0.03109070, 0.21370, 7.59570, 3.5876, 1.63820, 0.49294, 1},                     \
   {0.01554535, 0.20548, 14.1189, 6.1977, 3.36620, 0.62517, 1},                     \
   {0.01688690, 0.11125, 10.3570, 3.6231, 0.88026, 0.49671, 1}};

// This is the pw92 epsilon using the most accuracte parameters,
// and exact values for the exact constants.
template <typename num> static num pw92eps(const densvars<num> & d) {
  const parameter TUVWXYP[3][7] = PW92C_PARAMS;
#ifndef XCFUN_REF_PW92C
  const parameter c = 8.0 / (9.0 * (2 * pow(2, 1.0 / 3.0) - 2));
#else
  const parameter c = 1.709921;
#endif
  num zeta4 = pow(d.zeta, 4);
  num omegaval = omega(d.zeta);
  num sqrtr = sqrt(d.r_s);
  num e0 = eopt(sqrtr, TUVWXYP[0]);
  return e0 - eopt(sqrtr, TUVWXYP[2]) * omegaval * (1 - zeta4) / c +
         (eopt(sqrtr, TUVWXYP[1]) - e0) * omegaval * zeta4;
}

template <typename num> static num pw92eps_polarized(const num & a) {
  const parameter TUVWXYP[3][7] = PW92C_PARAMS;
  num sqrt_r_s = pow(3 / (4 * M_PI * a), 1.0 / 6.0);
  return eopt(sqrt_r_s, TUVWXYP[1]);
}
} // namespace pw92eps
