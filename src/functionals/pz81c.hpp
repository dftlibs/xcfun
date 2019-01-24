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
#include "functional.hpp"

namespace pz81eps {

// (1+zeta)^4/3 + (1-zeta)^(4/3) reorganized, plus constants.
template <typename num> static num fz(const densvars<num> & d) {
  const parameter p = pow(2.0, 4.0 / 3.0);
  const parameter q = (2 * pow(2.0, 1.0 / 3.0) - 2);
  return (p * (d.a_43 + d.b_43) * d.n_m13 / d.n - 2) / q;
}

template <typename num> static num Eld(const num & x, const parameter CB1B2[]) {
  return CB1B2[0] / (1 + CB1B2[1] * sqrt(x) + CB1B2[2] * x);
}

template <typename num> static num Ehd(const num & x, const parameter c[]) {
  return c[1] + log(x) * (c[0] + x * c[2]) + c[3] * x;
}

template <typename num> static num pz81eps(const densvars<num> & d) {
  parameter c[4][4] = {{-0.1423, 1.0529, 0.3334, 0},
                       {-0.0843, 1.3981, 0.2611, 0},
                       {0.0311, -0.048, 0.0020, -0.0116},
                       {0.01555000000, -0.0269, 0.0007, -0.0048}};
  if (1 > d.r_s)
    return Ehd(d.r_s, c[2]) + (Ehd(d.r_s, c[3]) - Ehd(d.r_s, c[2])) * fz(d);
  else
    return Eld(d.r_s, c[0]) + (Eld(d.r_s, c[1]) - Eld(d.r_s, c[0])) * fz(d);
}
}; // namespace pz81eps
