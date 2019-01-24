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
#include "slater.hpp"

FUNCTIONAL(XC_SLATERX) = {
    "Slater LDA exchange",
    "LDA Exchange functional\n"
    "P.A.M. Dirac, Proceedings of the Cambridge Philosophical "
    "Society, 26 (1930) 376.\n"
    "F. Bloch, Zeitschrift fuer Physik, 57 (1929) 545.\n\n"
    "Implemented by Ulf Ekstrom\n"
    "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
    XC_DENSITY,
    ENERGY_FUNCTION(slaterx) XC_A_B,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-11,
    {0.39E+02, 0.38E+02},
    {-0.241948147838E+03, // energy
     -0.420747936684E+01, // gradient
     -0.417120618800E+01,
     -0.359613621097E-01, // hessian
     0,
     -0.365895279649E-01}};
