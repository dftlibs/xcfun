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
#include "tpssc_eps.hpp"

template <typename num> static num tpssc(const densvars<num> & d) {
  num eps = tpssc_eps::tpssc_eps(d);
  return d.n * eps;
}

FUNCTIONAL(XC_TPSSC) = {"TPSS original correlation functional",
                        "TPSS original correlation functional.\n"
                        "J. Tao, J.P. Perdew, V. N. Staroverov, G. E. Scuseria,\n"
                        "Climbing the Density Functional Ladder:\n"
                        "Nonempirical Meta-Generalized Gradient Approximation\n"
                        "Designed for Molecules and Solids,\n"
                        "Phys. Rev. Lett. 91 (2003) 146401\n"
                        "Implemented by Andrea Debnarova\n",
                        XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                        ENERGY_FUNCTION(tpssc) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                        XC_PARTIAL_DERIVATIVES,
                        1,
                        1e-6,
                        {1, 2, 3, 4, 5, 6, 7},
                        {-2.1824017471364521e-01,
                         -0.114815778042036,
                         -7.968561473875205E-002,
                         8.304923723228601E-004,
                         1.671559231642408E-003,
                         8.327578102972385E-004,
                         -1.440645713140351E-005,
                         -1.440645713140351E-005}};

/* // Test case from A**
const double d[] = {0.153652558932587,
                    0.153652558932587,
                    8.390981882024583E-002,
                    8.390981882024583E-002,
                    8.390981882024583E-002,
                    6.826262722466284E-002,
                    6.826262722466284E-002};
const double ref[] = {0,//Fill in energy
                         -7.299596527683171E-002,
                     -7.299596527683171E-002,
                     1.247652223035057E-002, // Probably incorrect
                     5.978495466748646E-002, // Probably incorrect
                     1.247652223035057E-002, // Probably incorrect
                     -4.068185412322921E-002,
                     -4.068185412322921E-002};
*/
