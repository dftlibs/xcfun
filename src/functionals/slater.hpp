/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2018 Ulf Ekström and contributors.
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

#ifndef SLATER_H
#define SLATER_H

#include "constants.hpp"

template<class num>
static num slaterx(const densvars<num> &d) 
{ 
  return (-xc_constants::c_slater)*(d.a_43 + d.b_43);
}

#endif
