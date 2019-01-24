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

// M05 exchange functional. to be used with HF exchange factor of .28

template <typename num> static num m05x(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::fw;

  const parameter param_a[12] = {1.000000e+00,
                                 8.151000e-02,
                                 -4.395600e-01,
                                 -3.224220e+00,
                                 2.018190e+00,
                                 8.794310e+00,
                                 -2.950000e-03,
                                 9.820290e+00,
                                 -4.823510e+00,
                                 -4.817574e+01,
                                 3.648020e+00,
                                 3.402248e+01};

  return (pbex::energy_pbe_ab(pbex::R_pbe, d.a, d.gaa) * fw(param_a, d.a, d.taua) +
          pbex::energy_pbe_ab(pbex::R_pbe, d.b, d.gbb) * fw(param_a, d.b, d.taub));
}

FUNCTIONAL(XC_M05X) = {"M05 exchange",
                       "M05 Meta-Hybrid Exchange Functional\n"
                       "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory "
                       "Comput. 2, 364 (2006)\n"
                       "Implemented by Andre Gomes\n",
                       XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                       ENERGY_FUNCTION(m05x) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
                       XC_PARTIAL_DERIVATIVES,
                       1,
                       3e-5,
                       {1., .8, 1., 1., 1., 0.165, 0.1050},
                       {-1.57876583,
                        -2.12127045,
                        -2.11264351,
                        -0.00315462,
                        0.00000000,
                        -0.00444560,
                        3.45640232,
                        4.4349756}};
