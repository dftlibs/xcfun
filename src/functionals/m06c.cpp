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

template <typename num> static num m06c(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::m06_c_anti;
  using m0xy_metagga_xc_internal::m06_c_para;
  using m0xy_metagga_xc_internal::ueg_c_anti;
  using m0xy_metagga_xc_internal::ueg_c_para;
  using m0xy_metagga_xc_internal::zet;
  using pw91_like_x_internal::chi2;

  // parameters for anti-parallel spin contributions
  const parameter param_c_anti[5] = {
      3.741539e+00, 2.187098e+02, -4.531252e+02, 2.936479e+02, -6.287470e+01};
  const parameter param_d_anti[6] = {-2.741539e+00,
                                     -6.720113e-01,
                                     -7.932688e-02,
                                     1.918681e-03,
                                     -2.032902e-03,
                                     0.000000e+00};

  // parameters for parallel spin contributions
  const parameter param_c_para[5] = {
      5.094055e-01, -1.491085e+00, 1.723922e+01, -3.859018e+01, 2.845044e+01};
  const parameter param_d_para[6] = {4.905945e-01,
                                     -1.437348e-01,
                                     2.357824e-01,
                                     1.871015e-03,
                                     -3.788963e-03,
                                     0.000000e+00};

  num chi_a2 = chi2(d.a, d.gaa);
  num chi_b2 = chi2(d.b, d.gbb);
  num zet_a = zet(d.a, d.taua);
  num zet_b = zet(d.b, d.taub);
  num Dsigma_a = m0xy_metagga_xc_internal::Dsigma(d.a, d.gaa, d.taua);
  num Dsigma_b = m0xy_metagga_xc_internal::Dsigma(d.b, d.gbb, d.taub);

  // About six correct digits in Ec_ab
  num Ec_ab = ueg_c_anti(d) *
              m06_c_anti(param_c_anti, param_d_anti, chi_a2, zet_a, chi_b2, zet_b);
  num Ec_aa = ueg_c_para(d.a) *
              m06_c_para(param_c_para, param_d_para, chi_a2, zet_a, Dsigma_a);
  num Ec_bb = ueg_c_para(d.b) *
              m06_c_para(param_c_para, param_d_para, chi_b2, zet_b, Dsigma_b);
  return Ec_ab + Ec_aa + Ec_bb;
}

FUNCTIONAL(XC_M06C) = {
    "M06 Correlation",
    "M06 Meta-Hybrid Correlation Functional\n"
    "Zhao, Truhlar; Theor Chem Account (2008) 120:215-241\n"
    "Implemented by Andre Gomes\n"
    "Reference data from ADF, except energy which is self-computed.\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
    ENERGY_FUNCTION(m06c) XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
    XC_PARTIAL_DERIVATIVES,
    1,
    5e-5,
    {1, 2, 2, 1, 3, 1, 2},
    {-1.3563164150574752e-01,
     -0.131522634780454,
     -1.095804911824488E-002,
     -6.743582614131827E-004,
     0.000000000000000E+000,
     -7.065762242385136E-004,
     1.312970361640244E-002,
     -2.732803203907383E-003}};

/* Test from the midwest: */
/*const double d[] =
  {1., .8, 1., 1., 1., 0.165,   0.1050};
const double ref[] =
  { -1.57876583, -2.12127045, -2.11264351, -0.00315462, -0.00444560,  0.00000000,
1.72820116,  2.21748787 };
*/
/*
const double d[] = {0.153652558932587,
                           0.153652558932587,
                           8.390981882024769E-002,
                           8.390981882024769E-002,
                           8.390981882024769E-002,
                           6.826262722466429E-002,
                           6.826262722466429E-002};
const double ref[] = {-1.5172686665986701e-02, // Self-computed
                             -6.500066856941683E-002,
                             -6.500066856941683E-002,
                             -3.829432160900903E-002,
                             0.000000000000000E+000,
                             -3.829432160900903E-002,
                             8.813615476426437E-002,
                             8.813615476426437E-002};
*/
