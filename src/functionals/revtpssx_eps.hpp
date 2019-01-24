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

namespace revtpssx_eps {

template <typename num> static num epsx_unif(const num & d_n) {
  return -3 * cbrt(3 * PI2 * d_n) / (4 * PI);
}

template <typename num>
static num x(const num & d_n, const num & d_gnn, const num & d_tau) {
  parameter kapa = 0.804;
  parameter mu = 0.14;
  parameter b = 0.40;
  parameter e = 2.1677;
  parameter c = 2.35204;
  num p = d_gnn / (4 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
  num tauw = d_gnn / (8.0 * d_n);
  num z = tauw / d_tau;
  num z2 = pow2(z);
  num tau_unif = 0.3 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 5.0 / 3.0);
  num alpha = (d_tau - tauw) / tau_unif;
  num q_b = 9 * (alpha - 1) / (20 * sqrt(1 + b * alpha * (alpha - 1))) + 2 * p / 3.0;
  num x_a = p * (10.0 / 81.0 + c * z2 * z / pow2(1 + z2));
  x_a += 146.0 * pow2(q_b) / 2025.0;
  // Take the gradient out of the square root
  x_a -= 73.0 / 405.0 * q_b * d_gnn *
         sqrt((0.5 * 0.6 * 0.6) * pow(8 * d_n * d_tau, -2) + 0.5 * p * p);
  // instead of x_a -= 73.0*q_b*sqrt(0.5*0.6*0.6*z2 + 0.5*p*p)/405.0;
  x_a += pow2(10.0 * p / 81.0) / kapa;
  x_a += sqrt(e) * 0.6 * 0.6 * z2 * 20.0 / 81.0 + e * mu * pow3(p);
  return x_a / (1 + sqrt(e) * p);
}

template <typename num>
static num F_x(const num & d_n, const num & d_gnn, const num & d_tau) {
  parameter kapa = 0.804;
  num xpz = x(d_n, d_gnn, d_tau);
  return 1 + kapa - kapa / (1 + xpz / kapa);
}

template <typename num> static num revtpssx_eps(const densvars<num> & d) {
  num Fxa = F_x(2 * d.a, 4 * d.gaa, 2 * d.taua);
  num epsxunif_a = epsx_unif(2 * d.a);
  num Fxb = F_x(2 * d.b, 4 * d.gbb, 2 * d.taub);
  num epsxunif_b = epsx_unif(2 * d.b);
  return (epsxunif_a * Fxa * d.a + epsxunif_b * Fxb * d.b);
}
} // namespace revtpssx_eps
