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

#include "config.hpp"
#include "ctaylor.hpp"
#include "specmath.hpp"
#include "xcint.hpp"

#define FUNCTIONAL(F)                                                               \
  template <> const char * fundat_db<F>::symbol = #F;                               \
  template <> functional_data fundat_db<F>::d
#define EN(N, FUN) FUN<ctaylor<ireal_t, N>>,
#define ENERGY_FUNCTION(FUN) FOR_EACH(XC_MAX_ORDER, EN, FUN)
#define PARAMETER(P)                                                                \
  template <> const char * pardat_db<P>::symbol = #P;                               \
  template <> parameter_data pardat_db<P>::d
