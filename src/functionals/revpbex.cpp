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

#include "functional.hpp"
#include "pbex.hpp"

template <typename num> static num revpbex(const densvars<num> & d) {
  return pbex::energy_pbe_ab(pbex::R_revpbe, d.a, d.gaa) +
         pbex::energy_pbe_ab(pbex::R_revpbe, d.b, d.gbb);
}

FUNCTIONAL(XC_REVPBEX) = {"Revised PBE Exchange Functional",
                          "Revised PBE Exchange Functional\n"
                          "Y. Zhang and W., Phys. Rev. Lett 80, 890 (1998)\n"
                          "Implemented by Ulf Ekstrom and Andre Gomes\n",
                          XC_DENSITY | XC_GRADIENT,
                          ENERGY_FUNCTION(revpbex)};
