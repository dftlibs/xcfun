/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2020 Ulf Ekström and contributors.
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

#include <cmath>

#include "config.hpp"

#ifndef PI
constexpr auto PI = M_PI;
#endif

#ifndef PI2
constexpr auto PI2 = (M_PI * M_PI);
#endif

namespace xcfun_constants {
const parameter c_slater = pow(81 / (32 * M_PI), 1.0 / 3.0); // Typically called C_x
const parameter CF = 0.3 * pow(3 * PI2, 2.0 / 3.0);

// PBE constants.
const parameter param_gamma = (1 - log(2.0)) / (PI2);
constexpr parameter param_beta_pbe_paper = 0.066725;
constexpr parameter param_beta_accurate = 0.06672455060314922;
const parameter param_beta_gamma = param_beta_accurate / param_gamma;
} // namespace xcfun_constants
