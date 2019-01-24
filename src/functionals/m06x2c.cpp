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

// M06-2X correlation functional

template <typename num> static num m06x2c(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::m06_c_anti;
  using m0xy_metagga_xc_internal::m06_c_para;
  using m0xy_metagga_xc_internal::ueg_c_anti;
  using m0xy_metagga_xc_internal::ueg_c_para;
  using m0xy_metagga_xc_internal::zet;
  using pw91_like_x_internal::chi2;

  // parameters for anti-parallel spin contributions
  const parameter param_c_anti[5] = {
      8.833596e-01, 3.357972e+01, -7.043548e+01, 4.978271e+01, -1.852891e+01};
  const parameter param_d_anti[6] = {1.166404e-01,
                                     -9.120847e-02,
                                     -6.726189e-02,
                                     6.720580e-05,
                                     8.448011e-04,
                                     0.000000e+00};

  // parameters for parallel spin contributions
  const parameter param_c_para[5] = {
      3.097855e-01, -5.528642e+00, 1.347420e+01, -3.213623e+01, 2.846742e+01};
  const parameter param_d_para[6] = {6.902145e-01,
                                     9.847204e-02,
                                     2.214797e-01,
                                     -1.968264e-03,
                                     -6.775479e-03,
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

FUNCTIONAL(XC_M06X2C) = {"M06-2X Correlation",
                         "M06-2X Meta-Hybrid Correlation Functional\n"
                         "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory "
                         "Comput. 2, 364 (2006)\n"
                         "Implemented by Andre Gomes\n",
                         XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                         ENERGY_FUNCTION(m06x2c)};

/*  const double d[] =
  {1., .8, 1., 1., 1., 0.165,   0.1050}; //.33, .21};
const double ref[] =
  { -1.57876583, -2.12127045, -2.11264351, -0.00315462, -0.00444560,  0.00000000,
1.72820116,  2.21748787 };
f.add_test(XC_VARS_AB,1,d,ref,1e-5); */
