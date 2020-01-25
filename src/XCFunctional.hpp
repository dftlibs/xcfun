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

#include <array>

#include "XCFun/xcfun.h"
#include "functionals/list_of_functionals.hpp"

struct functional_data;

struct XCFunctional {
  explicit XCFunctional();

  int nr_active_functionals{0};
  int order{-1};
  int depends{0}; // XC_DENSITY, gradient etc
  xcfun_mode mode{XC_MODE_UNSET};
  xcfun_vars vars{XC_VARS_UNSET};
  std::array<functional_data *, XC_NR_FUNCTIONALS> active_functionals{nullptr};
  std::array<double, XC_NR_PARAMETERS_AND_FUNCTIONALS> settings;
};
