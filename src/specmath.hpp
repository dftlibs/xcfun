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

#ifdef XCFUN_NO_ERF
void xcint_die(const char * message, int code);
template <typename T> T erf(T x) {
  xcint_die("XcFun erf called but XCFUN_NO_ERF was defined", 0);
  return 0;
}
#endif

// Some math-related functions useful in many places
#include "ctaylor.hpp"

template <typename T> static T pow2(const T & t) { return t * t; }

template <typename T> static T pow3(const T & t) { return t * t * t; }

template <typename T, class T2>
static T poly(const T & x, int ndeg, const T2 coeffs[]) {
  // Horner rule
  T res = coeffs[--ndeg];
  while (ndeg) {
    res *= x;
    res += coeffs[--ndeg];
  }
  return res;
}

template <typename T, class S> static T ufunc(const T & x, S a) {
  return pow(1 + x, a) + pow(1 - x, a);
}
