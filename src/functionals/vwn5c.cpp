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
#include "vwn.hpp"

template <typename num> static num vwn5c(const densvars<num> & d) {
  return d.n * vwn::vwn5_eps(d);
}

FUNCTIONAL(XC_VWN5C) = {
    "VWN5 LDA Correlation functional",
    "VWN5 LDA Correlation functional\n"
    "S.H. Vosko, L. Wilk, and M. Nusair: Accurate spin-dependent\n"
    "electron liquid correlation energies for local spin density\n"
    "calculations: a critical analysis, Can. J. Phys. 58 (1980) 1200-1211.\n"
    "Originally from Dalton, polished and converted by Ulf Ekstrom.\n"
    "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_vwn5.html\n",
    XC_DENSITY,
    ENERGY_FUNCTION(vwn5c) XC_A_B,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-11,
    {0.39E+02, 0.38E+02},
    {-0.851077910672E+01,
     -0.119099058995E+00,
     -0.120906044904E+00,
     0.756836181702E-03,
     -0.102861281830E-02,
     0.800136175083E-03}};
