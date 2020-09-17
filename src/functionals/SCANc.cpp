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

#include "SCAN_like_eps.hpp"
#include "constants.hpp"
#include "functional.hpp"

template <class num> static num SCANC(const densvars<num> & d) {
  num eps_c = SCAN_eps::SCAN_C(d, 0, 0, 0);

  return eps_c;
}

FUNCTIONAL(XC_SCANC) = {
    "SCAN correlation functional",
    "SCAN correlation functional.\n"
    "J. Sun, A. Ruzsinszky, and J. P. Perdew, Phys. Rev. Lett. 115, 036402 (2015)."
    "Implemented by James Furness (@JFurness1)\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
    ENERGY_FUNCTION(SCANC) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
    XC_PARTIAL_DERIVATIVES,
    1,
    1e-11,
    {0.217, 0.0632, 0.191, 0.0535, 0.015, 0.267, 0.0328},
    {-7.42689151900e-03,
     -1.91882337268e-02,
     -6.34666733088e-02,
     1.08870558794e-02,
     2.17741117588e-02,
     1.08870558794e-02,
     -1.76437503720e-02,
     -1.76437503720e-02}};
