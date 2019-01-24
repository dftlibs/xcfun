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

// Use #define XC_NO_REGULARIZATION to turn off
// checks and balances for physical densities.

// Enable functionals still in development (probably buggy)
#define XCFUN_IN_DEVELOPMENT

// Use #define XCFUN_REF_PW92C to use inaccurate constants in
// PW92C. This matches the reference implementation.

// Use inaccurate mu value in pbe exchange.
// #define XCFUN_REF_PBEX_MU

// This is the internal scalar type of the library, can be
// different from the external interface.
#ifndef WITH_QD
typedef double ireal_t;
#define INNER_TO_OUTER(INNER) INNER
#else
#include <qd/qd_real.h>
typedef qd_real ireal_t;
#define XCFUN_NUM_CONVERT // Must convert real types at i/o
#define INNER_TO_OUTER(INNER) to_double(INNER)
#endif
