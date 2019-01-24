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
#include "m0xy_fun.hpp"
#include "pbex.hpp"

// M05-2X exchange functional. to be used with HF exchange factor of .56

template <typename num> static num m05x2x(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::fw;

  const parameter param_a[12] = {1.000000e+00,
                                 -5.683300e-01,
                                 -1.300570e+00,
                                 5.500700e+00,
                                 9.064020e+00,
                                 -3.221075e+01,
                                 -2.373298e+01,
                                 7.022996e+01,
                                 2.988614e+01,
                                 -6.025778e+01,
                                 -1.322205e+01,
                                 1.523694e+01};

  return (pbex::energy_pbe_ab(pbex::R_pbe, d.a, d.gaa) * fw(param_a, d.a, d.taua) +
          pbex::energy_pbe_ab(pbex::R_pbe, d.b, d.gbb) * fw(param_a, d.b, d.taub));
}

FUNCTIONAL(XC_M05X2X) = {"M05-2X exchange",
                         "M05-2X Meta-Hybrid Exchange Functional\n"
                         "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory "
                         "Comput. 2, 364 (2006)\n"
                         "Implemented by Andre Gomes\n",
                         XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                         ENERGY_FUNCTION(m05x2x) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                         XC_PARTIAL_DERIVATIVES,
                         1,
                         3e-5,
                         {1., .8, 1., 1., 1., 0.165, 0.1050},
                         {-1.38233309,
                          -0.19638222,
                          -0.08614105,
                          -0.00289174,
                          0.00000000,
                          -0.00365982,
                          -3.18842316000000,
                          -3.90587738}};
