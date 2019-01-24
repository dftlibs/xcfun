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

// M05 correlation functional

template <typename num> static num m05c(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::m05_c_anti;
  using m0xy_metagga_xc_internal::m05_c_para;
  using m0xy_metagga_xc_internal::ueg_c_anti;
  using m0xy_metagga_xc_internal::ueg_c_para;
  using m0xy_metagga_xc_internal::zet;
  using pw91_like_x_internal::chi2;

  // parameters for anti-parallel spin contributions
  const parameter param_c_anti[5] = {
      1.000000e+00, 3.785690e+00, -1.415261e+01, -7.465890e+00, 1.794491e+01};

  // parameters for parallel spin contributions
  const parameter param_c_para[5] = {
      1.000000e+00, 3.773440e+00, -2.604463e+01, 3.069913e+01, -9.226950e+00};

  num chi_a2 = chi2(d.a, d.gaa);
  num chi_b2 = chi2(d.b, d.gbb);
  num zet_a = zet(d.a, d.taua);
  num zet_b = zet(d.b, d.taub);
  num Dsigma_a = m0xy_metagga_xc_internal::Dsigma(d.a, d.gaa, d.taua);
  num Dsigma_b = m0xy_metagga_xc_internal::Dsigma(d.b, d.gbb, d.taub);

  num Ec_ab = ueg_c_anti(d) * m05_c_anti(param_c_anti, chi_a2, chi_b2);
  num Ec_aa = ueg_c_para(d.a) * m05_c_para(param_c_para, chi_a2, zet_a, Dsigma_a);
  num Ec_bb = ueg_c_para(d.b) * m05_c_para(param_c_para, chi_b2, zet_b, Dsigma_b);

  return Ec_ab + Ec_aa + Ec_bb;
}

FUNCTIONAL(XC_M05C) = {"M05 Correlation",
                       "M05 Meta-Hybrid Correlation Functional\n"
                       "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory "
                       "Comput. 2, 364 (2006)\n"
                       "Implemented by Andre Gomes\n",
                       XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                       ENERGY_FUNCTION(m05c)};
/*  const double d[] =
  {1., .8, 1., 1., 1., .33, .21};
const double ref[] =
  { -0.06599246, -0.15418438,  0.02729798,  0.03090769,  0.04788146,  0.00000000,
-0.09757618, -0.23742358 };
f.add_test(XC_VARS_AB,1,d,ref,2e-5); */
