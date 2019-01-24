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

#include "b97xc.hpp"
#define PREFACTOR 0.9305257363491002 // 1.5*pow(3/(4*PI),1.0/3.0)

namespace b97x {

// LSDA factor, for alpha and beta spin

template <typename num> static num e_x_LSDA_ab(const num & a_43) {
  return -PREFACTOR * a_43;
}

// parameters for enhancement factor; c0,c1,c2

const parameter c_b97[3] = {0.8094, 0.5073, 0.7481};
const parameter c_b97_1[3] = {0.789518, 0.573805, 0.660975};
const parameter c_b97_2[3] = {0.827642, 0.047840, 1.76125};
const parameter Gamma = 0.004;

template <typename num>
static num energy_b97x_ab(const parameter & Gamma,
                          const parameter c_params[],
                          const num & a_43,
                          const num & gaa) {
  num s2_ab = b97xc::spin_dens_gradient_ab2(gaa, a_43);
  return e_x_LSDA_ab(a_43) * b97xc::enhancement(Gamma, c_params, s2_ab);
}
} // namespace b97x
