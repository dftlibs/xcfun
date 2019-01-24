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
#include "pbec_eps.hpp"

namespace tpssx_eps {

template <typename num> static num fx_unif(const num & d) {
  return (-0.75 * pow(3 / PI, 1.0 / 3.0)) * pow(d, 4.0 / 3.0);
}

template <typename num>
static num x(const num & d_n, const num & d_gnn, const num & d_tau) {
  const parameter kappa = 0.804;
  const parameter mu = 0.21951;
  const parameter b = 0.40;
  const parameter e = 1.537;
  const parameter c = 1.59096;
  num p0 = 1.0 / (4 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
  num p = d_gnn * p0;
  num tauw = d_gnn / (8.0 * d_n);
  num z = tauw / d_tau;
  num z2 = pow2(z);
  num tau_unif = 0.3 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 5.0 / 3.0);
  num alpha = (d_tau - tauw) / tau_unif; //  5*p*z/(3.0 - 3.0*z);
  num q_b =
      (9.0 / 20.0) * (alpha - 1) / sqrt(1 + b * alpha * (alpha - 1)) + 2 * p / 3.0;
  num x_a = p * (10.0 / 81.0 + c * z2 / pow2(1 + z2));
  x_a += 146.0 * pow2(q_b) / 2025.0;
  //    x_a -= 73.0/405.0*q_b*sqrt(0.5*0.6*0.6*z2 + 0.5*p*p);
  // Take the gradient out of the square root
  x_a -= 73.0 / 405.0 * q_b * d_gnn *
         sqrt((0.5 * 0.6 * 0.6) * pow(8 * d_n * d_tau, -2) + 0.5 * p0 * p0);
  x_a += pow2(p * 10.0 / 81.0) / kappa;
  x_a += 2 * sqrt(e) * 0.6 * 0.6 * z2 * 10.0 / 81.0 + e * mu * pow3(p);
  return x_a / pow2(1 + sqrt(e) * p);
}

template <typename num>
static num F_x(const num & d_n, const num & d_gnn, const num & d_tau) {
  const parameter kappa = 0.804;
  num xpz = x(d_n, d_gnn, d_tau);
  return 1 + kappa - kappa / (1 + xpz / kappa);
}
} // namespace tpssx_eps
