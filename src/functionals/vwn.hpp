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

#include "functional.hpp"

namespace vwn {

static inline parameter vwn_a(const parameter p[]) {
  return p[0] * p[2] / (p[0] * p[0] + p[0] * p[2] + p[3]) - 1;
}

static inline parameter vwn_b(const parameter p[]) {
  return 2 * (p[0] * p[2] / (p[0] * p[0] + p[0] * p[2] + p[3]) - 1) + 2;
}

static inline parameter vwn_c(const parameter p[]) {
  return 2 * p[2] *
         (1 / sqrt(4 * p[3] - p[2] * p[2]) -
          p[0] / ((p[0] * p[0] + p[0] * p[2] + p[3]) * sqrt(4 * p[3] - p[2] * p[2]) /
                  (p[2] + 2 * p[0])));
}

template <typename num> static num vwn_x(const num & s, const parameter p[]) {
  return s * s + p[2] * s + p[3];
}

template <typename num> static num vwn_y(const num & s, const parameter p[]) {
  return s - p[0];
}

template <typename num> static num vwn_z(const num & s, const parameter p[]) {
  return sqrt(4 * p[3] - p[2] * p[2]) / (2 * s + p[2]);
}

template <typename num> static num vwn_f(const num & s, const parameter p[]) {
  return 0.5 * p[1] *
         (2 * log(s) + vwn_a(p) * log(vwn_x(s, p)) - vwn_b(p) * log(vwn_y(s, p)) +
          vwn_c(p) * atan(vwn_z(s, p)));
}

template <typename num> static num vwn5_eps(const densvars<num> & d) {
// ulfek: second elements are multiplied by 2 wrt molpro manual
#ifndef XCFUN_VWN5_REF
  const parameter para[] = {-0.10498, 0.0621814, 3.72744, 12.9352};
  const parameter ferro[] = {-0.325, 0.0310907, 7.06042, 18.0578};
  const parameter inter[] = {
      -0.0047584, -pow(3 * M_PI * M_PI, -1.0), 1.13107, 13.0045};
#else
  const parameter para[] = {
      -0.10498, 0.0621813817393097900698817274255, 3.72744, 12.9352};
  const parameter ferro[] = {
      -0.325, 0.0310906908696548950349408637127, 7.06042, 18.0578};
  const parameter inter[] = {
      -0.0047584, -pow(3 * M_PI * M_PI, -1.0), 1.13107, 13.0045};
#endif
  num s = sqrt(d.r_s);
  // Constant is (2^1/3-1)^-1/2
  num g = 1.92366105093154 * (ufunc(d.zeta, 4.0 / 3.0) - 2);
  num zeta4 = pow(d.zeta, 4);
  // FIXME: 1 - zeta^4 has a cancellation free form (ask maxima)
  num dd =
      g * ((vwn_f(s, ferro) - vwn_f(s, para)) * zeta4 +
           vwn_f(s, inter) * (1 - zeta4) * (9.0 / 4.0 * (pow(2, 1.0 / 3.0) - 1)));
  return (vwn_f(s, para) + dd);
}

template <typename num> static num vwn3_eps(const densvars<num> & d) {
  // ulfek: This table taken verbatim from Dalton
  const parameter para[] = {-0.4092860, 0.0621814, 13.0720, 42.7198};
  const parameter ferro[] = {-0.7432940, 0.0310907, 20.1231, 101.578};
  const parameter inter[] = {-0.0047584, -0.0337737, 1.13107, 13.0045};
  num s = sqrt(d.r_s);
  // Constant is (2^1/3-1)^-1/2
  num g = 1.92366105093154 * (ufunc(d.zeta, 4.0 / 3.0) - 2);
  num zeta4 = pow(d.zeta, 4);
  // FIXME: 1 - zeta^4 has a cancellation free form (ask maxima)
  num dd =
      g * ((vwn_f(s, ferro) - vwn_f(s, para)) * zeta4 +
           vwn_f(s, inter) * (1 - zeta4) * (9.0 / 4.0 * (pow(2, 1.0 / 3.0) - 1)));
  return (vwn_f(s, para) + dd);
}
} // namespace vwn
