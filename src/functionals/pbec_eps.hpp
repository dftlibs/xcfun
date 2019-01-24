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
#include "pw92eps.hpp"
#include "vwn.hpp"

namespace pbec_eps {
template <typename num, class T> static num A(const num & eps, const T & u3) {
  using xc_constants::param_beta_gamma;
  using xc_constants::param_gamma;
  return param_beta_gamma / expm1(-eps / (param_gamma * u3));
}

template <typename num, class T>
static num H(const num & d2, const num & eps, const T & u3) {
  num d2A = d2 * A(eps, u3);
  using xc_constants::param_beta_gamma;
  using xc_constants::param_gamma;
  return param_gamma * u3 *
         log(1 + param_beta_gamma * d2 * (1 + d2A) / (1 + d2A * (1 + d2A)));
}

// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template <typename num> static num phi(const densvars<num> & d) {
  return pow(2.0, -1.0 / 3.0) * d.n_m13 * d.n_m13 * (sqrt(d.a_43) + sqrt(d.b_43));
}

template <typename num> static num pbec_eps(const densvars<num> & d) {
  num eps = pw92eps::pw92eps(d);
  num u = phi(d);
  // Avoiding the square root of d.gnn here
  num d2 = pow(1.0 / 12 * pow(3, 5.0 / 6.0) / pow(M_PI, -1.0 / 6), 2.0) * d.gnn /
           (u * u * pow(d.n, 7.0 / 3.0));
  return (eps + H(d2, eps, pow3(u)));
}

template <typename num>
static num pbec_eps_polarized(const num & a, const num & gaa) {
  num eps = pw92eps::pw92eps_polarized(a);
  parameter u = pow(2.0, -1.0 / 3.0); // phi(d) for alpha or beta density =0

  // Avoiding the square root of d.gnn here
  num d2 = pow(1.0 / 12 * pow(3, 5.0 / 6.0) / pow(M_PI, -1.0 / 6), 2.0) * gaa /
           (u * u * pow(a, 7.0 / 3.0));
  return (eps + H(d2, eps, pow3(u)));
}
} // namespace pbec_eps
