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

// M06-2X exchange functional. to be used with HF exchange factor of .54
// because the param_d[5] array is all zero, it is not included here, and
// therefore the lsda_x() * h() terms (see e.g M06x) drop out, as h()=0

template <typename num> static num m06x2x(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::fw;

  const parameter param_a[12] = {4.600000e-01,
                                 -2.206052e-01,
                                 -9.431788e-02,
                                 2.164494e+00,
                                 -2.556466e+00,
                                 -1.422133e+01,
                                 1.555044e+01,
                                 3.598078e+01,
                                 -2.722754e+01,
                                 -3.924093e+01,
                                 1.522808e+01,
                                 1.522227e+01};

  return (pbex::energy_pbe_ab(pbex::R_pbe, d.a, d.gaa) * fw(param_a, d.a, d.taua) +
          pbex::energy_pbe_ab(pbex::R_pbe, d.b, d.gbb) * fw(param_a, d.b, d.taub));
}

FUNCTIONAL(XC_M06X2X) = {
    "M06-2X exchange",
    "M06-2X Meta-Hybrid Exchange Functional\n"
    "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
    "Implemented by Andre Gomes\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
    ENERGY_FUNCTION(m06x2x) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
    XC_PARTIAL_DERIVATIVES,
    1,
    3e-5,
    {1., .8, 1., 1., 1., 0.165, 0.1050},
    {-0.63803890,
     -0.81863653,
     -0.81208750,
     -0.00127795,
     0.00000000,
     -0.00179117,
     1.25220996,
     1.60808316}};
