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

// M06-L exchange functional. does not use HF exchange

template <typename num> static num m06lx(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::alpha_x;
  using m0xy_metagga_xc_internal::fw;
  using m0xy_metagga_xc_internal::h;
  using m0xy_metagga_xc_internal::zet;
  using pw91_like_x_internal::chi2;

  const parameter param_a[12] = {3.987756e-01,
                                 2.548219e-01,
                                 3.923994e-01,
                                 -2.103655e+00,
                                 -6.302147e+00,
                                 1.097615e+01,
                                 3.097273e+01,
                                 -2.318489e+01,
                                 -5.673480e+01,
                                 2.160364e+01,
                                 3.421814e+01,
                                 -9.049762e+00};
  const parameter param_d[6] = {6.012244e-01,
                                4.748822e-03,
                                -8.635108e-03,
                                -9.308062e-06,
                                4.482811e-05,
                                0.000000e+00};

  num chia2 = chi2(d.a, d.gaa);
  num chib2 = chi2(d.b, d.gbb);

  return ((pbex::energy_pbe_ab(pbex::R_pbe, d.a, d.gaa) * fw(param_a, d.a, d.taua) +
           lsda_x(d.a) * h(param_d, alpha_x, chia2, zet(d.a, d.taua))) +
          (pbex::energy_pbe_ab(pbex::R_pbe, d.b, d.gbb) * fw(param_a, d.b, d.taub) +
           lsda_x(d.b) * h(param_d, alpha_x, chib2, zet(d.b, d.taub))));
}

FUNCTIONAL(XC_M06LX) = {
    "M06-L exchange",
    "M06-L Meta-GGA Exchange Functional\n"
    "Y Zhao and D. G. Truhlar, J. Chem. Phys. 125, 194101 (2006)\n"
    "Implemented by Andre Gomes\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
    ENERGY_FUNCTION(m06lx) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
    XC_PARTIAL_DERIVATIVES,
    1,
    1e-5,
    {1., .8, 1., 1., 1., 0.165, 0.1050},
    {-1.60059999,
     -1.85109109,
     -1.81344820,
     -0.00370516,
     0.00000000,
     -0.00508827,
     2.3886147,
     3.05312252}};
