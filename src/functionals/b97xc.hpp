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

#include "xcint.hpp"

namespace b97xc {

// square of spin-density gradient
template <typename num>
static num spin_dens_gradient_ab2(const num & gaa, const num & a_43) {

  return abs(gaa) / a_43 / a_43;
}

template <typename num>
static num ux_ab(const parameter & Gamma, const num & spin_dens_grad) {

  return Gamma * spin_dens_grad / (1.0 + Gamma * spin_dens_grad);
}

template <typename num>
static num enhancement(const parameter & Gamma,
                       const parameter c_params[],
                       const num & spin_dens_grad) {
  num ux = ux_ab(Gamma, spin_dens_grad);

  return c_params[0] + c_params[1] * ux + c_params[2] * ux * ux;
}
} // namespace b97xc
