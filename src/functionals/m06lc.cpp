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

// M06 correlation functional

template <typename num> static num m06lc(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::m06_c_anti;
  using m0xy_metagga_xc_internal::m06_c_para;
  using m0xy_metagga_xc_internal::ueg_c_anti;
  using m0xy_metagga_xc_internal::ueg_c_para;
  using m0xy_metagga_xc_internal::zet;
  using pw91_like_x_internal::chi2;

  // parameters for anti-parallel spin contributions
  const parameter param_c_anti[5] = {
      6.042374e-1, 1.776783e2, -2.513252e2, 7.635173e1, -1.255699e1};
  const parameter param_d_anti[6] = {3.957626e-1,
                                     -5.614546e-1,
                                     1.403963e-2,
                                     9.831442e-4,
                                     -3.577176e-3,
                                     0.000000e+00};

  // parameters for parallel spin contributions
  const parameter param_c_para[5] = {
      5.349466e-1, 5.396620e-1, -3.161217e1, 5.149592e1, -2.919613e1};
  const parameter param_d_para[6] = {4.650534e-1,
                                     1.617589e-1,
                                     1.833657e-1,
                                     4.692100e-4,
                                     -4.990573e-3,
                                     0.000000e+00};

  num chi_a2 = chi2(d.a, d.gaa);
  num chi_b2 = chi2(d.b, d.gbb);
  num zet_a = zet(d.a, d.taua);
  num zet_b = zet(d.b, d.taub);
  num Dsigma_a = m0xy_metagga_xc_internal::Dsigma(d.a, d.gaa, d.taua);
  num Dsigma_b = m0xy_metagga_xc_internal::Dsigma(d.b, d.gbb, d.taub);

  num Ec_ab = ueg_c_anti(d) *
              m06_c_anti(param_c_anti, param_d_anti, chi_a2, zet_a, chi_b2, zet_b);
  num Ec_aa = ueg_c_para(d.a) *
              m06_c_para(param_c_para, param_d_para, chi_a2, zet_a, Dsigma_a);
  num Ec_bb = ueg_c_para(d.b) *
              m06_c_para(param_c_para, param_d_para, chi_b2, zet_b, Dsigma_b);
  return Ec_ab + Ec_aa + Ec_bb;
}

FUNCTIONAL(XC_M06LC) = {"M06-L Correlation",
                        "M06-L Meta GGA Correlation Functional\n"
                        "Zhao, Truhlar, JCP 125, 194101 (2006)\n"
                        "Implemented by Andre Gomes & Ulf Ekstrom\n",
                        XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                        ENERGY_FUNCTION(m06lc)};
