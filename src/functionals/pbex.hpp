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

#include "pw9xx.hpp"

namespace pbex {
// values of R for PBEx and REVPBEx

const parameter R_pbe = 0.804;
const parameter R_revpbe = 1.245;

// enhancement factor F(S), common to PBEx and REVPBEx
template <typename num>
static num enhancement(const parameter & R, const num & rho, const num & grad) {
  using pw91_like_x_internal::S2;

#ifdef XCFUN_REF_PBEX_MU
  // ulfek: mu from Daresbury implementation
  const parameter mu = 0.2195149727645171;
#else
  const parameter mu = 0.066725 * M_PI * M_PI / 3.0;
#endif
  num st2 = S2(rho, grad);
  num t1 = 1 + mu * st2 / R;
  return 1 + R - R / t1; // Intel <= 11.1 miscompiles(?) this line with -fast
}

template <typename num>
static num enhancement_RPBE(const num & rho, const num & grad) {
  using pw91_like_x_internal::S2;
  const parameter mu = 0.2195149727645171;
  return 1 - R_pbe * expm1((-mu / R_pbe) * S2(rho, grad));
}

template <typename num>
static num energy_pbe_ab(const parameter & R, const num & rho, const num & grad) {
  using pw91_like_x_internal::prefactor;
  return prefactor(rho) * enhancement(R, rho, grad);
}
} // namespace pbex
