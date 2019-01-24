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

// common functions for exchange (and kinetic energy) functionals
//
// related to PW91x - like PW91x (and PW91k, which is a reparametrization
// of PW91x to describe T_s), PBEx, PBEREVx and so on

namespace pw91_like_x_internal {
// formulas for the auxiliary quantities below can be found, for
// instance, at
//
//    http://www.molpro.net/info/current/doc/manual/node734.html#dftfun:PW91X
//

/*
  The chi2 and S2 functions return the square of the typical chi and S
  functions. This avoids a square root of the gradient square norm,
  which makes the functions differentiable at grad = 0. This is
  physical, so the square root can always(?) be avoided. A special
  sqrt(x)asinh(sqrt(x)) function is implemented for use in this context.
 */
template <typename num> static num chi2(const num & rho, const num & grad) {
  return grad / pow(rho, 8.0 / 3.0);
}

template <typename num> static num S2(const num & rho, const num & grad) {
  return grad / pow(rho, 8.0 / 3.0) *
         pow(pow(6.0, 2.0 / 3.0) / (12 * pow(M_PI, 2.0 / 3.0)), 2.0);
}

// prefactor multiples the enhancement factor F(S), which is then different
// for the different functionals

template <typename num> static num prefactor(const num & rho) {
  // aspg: the 2^.333 factor here i can't see in the molpro formula, will have to
  // double-check this - but as it is it matches the results for the database
  // e.g
  //   http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_pw91.html
  //   http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_pbe.html
  // and for pbe, it is close to the dalton code..(trunncation of variables
  // there make things deviate?)
  // ue: More likely the parameters, Dalton has some strange choice of
  // precision in its constants.
  return -0.75 * pow(2.0, 1.0 / 3.0) * pow(3 * M_PI * M_PI, 1.0 / 3.0) *
         pow(rho, 4.0 / 3.0) / M_PI;
}

// prefactor for the pw91k functional
template <typename num> static num pw91k_prefactor(const num & rho) {
  using xc_constants::CF;

  return CF * pow(2.0, 2.0 / 3.0) * pow(rho, 5.0 / 3.0);
}

// enhancement factor F(S), common to PW91x and PW91k
template <typename num>
static num pw91xk_enhancement(const parameter param_AB[6],
                              const num & rho,
                              const num & grad) {
  // This formula contains a square root of grad, which is not differentiable at
  // zero.
  // num st = S(rho,grad);
  // num t1 = 1 + param_AB[0]*st*asinh(param_AB[1]*st);

  // The method below never takes the square root of S explicitly,
  // and can deal with zero gradients.

  num st2 = S2(rho, grad);
  num t1 =
      1 + param_AB[0] * sqrtx_asinh_sqrtx(pow(param_AB[1], 2) * st2) / param_AB[1];
  num t2 = st2 * (param_AB[2] - param_AB[3] * exp(-param_AB[4] * st2));

  num numerator = t1 + t2;
  num denominator = t1 + param_AB[5] * st2 * st2;

  return numerator / denominator;
}
} // namespace pw91_like_x_internal
