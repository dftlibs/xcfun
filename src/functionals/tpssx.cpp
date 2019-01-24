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

#include "constants.hpp"
#include "functional.hpp"
#include "tpssx_eps.hpp"

using namespace tpssx_eps;

template <typename num> static num tpssx(const densvars<num> & d) {
  num Fxa = F_x(2 * d.a, 4 * d.gaa, 2 * d.taua);
  num epsxunif_a = fx_unif(2 * d.a);
  num Fxb = F_x(2 * d.b, 4 * d.gbb, 2 * d.taub);
  num epsxunif_b = fx_unif(2 * d.b);
  return 0.5 * (epsxunif_a * Fxa + epsxunif_b * Fxb);
}

FUNCTIONAL(XC_TPSSX) = {"TPSS original exchange functional",
                        "TPSS original exchange functional.\n"
                        "J. Tao, J.P. Perdew, V. N. Staroverov, G. E. Scuseria,\n"
                        "Climbing the Density Functional Ladder:\n"
                        "Nonempirical Meta-Generalized Gradient Approximation\n"
                        "Designed for Molecules and Solids,\n"
                        "Phys. Rev. Lett. 91 (2003) 146401\n"
                        "Implemented by Andrea Debnarova\n",
                        XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                        ENERGY_FUNCTION(tpssx) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                        XC_PARTIAL_DERIVATIVES,
                        1,
                        1e-8,
                        {0.153652558932587,
                         0.153652558932587,
                         8.390981882024583E-002,
                         8.390981882024583E-002,
                         8.390981882024583E-002,
                         6.826262722466284E-002,
                         6.826262722466284E-002},
                        {-1.7369862616505458e-01, // XCFun energy from ADF
                         -0.696647334686932,
                         -0.696647334686932,
                         -0.104380091839662,
                         0.0,
                         -0.104380091839662,
                         0.128315600541666,
                         0.128315600541666}};
