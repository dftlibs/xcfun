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

#include <cstdio>
#include <cstdlib>

namespace xcfun {
/*! Used for regularizing input */
constexpr auto XCFUN_TINY_DENSITY = 1e-14;

inline void die(const char * message, int code) {
  std::fprintf(stderr, "XCFun fatal error %i: ", code);
  std::fprintf(stderr, "%s", message);
  std::fprintf(stderr, "\n");
  std::exit(-1);
}
} // namespace xcfun

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

typedef ireal_t parameter;
