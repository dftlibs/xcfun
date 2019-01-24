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

#ifndef M_PI // M_PI is not standard for some reason
#define M_PI 3.14159265358979323846
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef PI2
#define PI2 (M_PI * M_PI)
#endif

namespace xc_constants {
const parameter c_slater = pow(81 / (32 * M_PI), 1.0 / 3.0); // Typically called C_x
const parameter CF = 0.3 * pow(3 * M_PI * M_PI, 2.0 / 3.0);

// PBE constants.
const parameter param_gamma = (1 - log(2.0)) / (M_PI * M_PI);
const parameter param_beta_pbe_paper = 0.066725;
const parameter param_beta_accurate = 0.06672455060314922;
const parameter param_beta_gamma = param_beta_accurate / param_gamma;
} // namespace xc_constants
