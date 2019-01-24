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

template <typename num> static num becke_alpha(const num & na, const num & gaa) {
  const parameter c = pow(81 / (4 * M_PI), 1.0 / 3.0) / 2;
  const parameter d = 0.0042;
  num na43 = pow(na, 4.0 / 3.0);
  num lda = -c * na43;
  num chi2 = gaa * pow(na, -8.0 / 3.0);
  num b88 = -(d * na43 * chi2) / (1 + 6 * d * sqrtx_asinh_sqrtx(chi2));
  return lda + b88;
}

template <typename num> static num becke_corr(const num & na, const num & gaa) {
  const parameter d = 0.0042;
  num na43 = pow(na, 4.0 / 3.0);
  num chi2 = gaa * pow(na, -8.0 / 3.0);
  return -(d * na43 * chi2) / (1 + 6 * d * sqrtx_asinh_sqrtx(chi2));
}

// Short range becke exchange as used in camb3lyp If mu=0 this reduces
// to the standard beckex for which must be used instead.
// FIXME: The erf + something is basically the erf minus
// its asymptotic expansion. This is horribel for numerics,
// will have to code a special function.
// As coded here the code will fail if mu = 0, in which case the
// regular beckex should be used.

template <typename num>
static num becke_sr(parameter mu, const num & na, const num & gaa) {
  const parameter cparam = pow(81 / (4 * M_PI), 1.0 / 3.0) / 2;
  const parameter d = 0.0042;
  num na43 = pow(na, 4.0 / 3.0);
  num chi2 = gaa * pow(na, -8.0 / 3.0);
  num K = 2 * (cparam + (d * chi2) / (1 + 6 * d * sqrtx_asinh_sqrtx(chi2)));
  num a = mu * sqrt(K) / (6 * sqrt(M_PI) * pow(na, 1.0 / 3.0));
  num b = expm1(-1 / (4 * a * a));
  num c = 2 * a * a * b + 0.5;
  return -0.5 * na43 * K *
         (1 - 8.0 / 3.0 * a * (sqrt(M_PI) * erf(1 / (2 * a)) + 2 * a * (b - c)));
}

template <typename num>
static num becke_cam(parameter alpha,
                     parameter beta,
                     parameter mu,
                     const num & na,
                     const num & gaa) {
  const parameter cparam = pow(81 / (4 * M_PI), 1.0 / 3.0) / 2;
  const parameter d = 0.0042;
  num na43 = pow(na, 4.0 / 3.0);
  num chi2 = gaa * pow(na, -8.0 / 3.0);
  num K = 2 * (cparam + (d * chi2) / (1 + 6 * d * sqrtx_asinh_sqrtx(chi2)));
  num a = mu * sqrt(K) / (6 * sqrt(M_PI) * pow(na, 1.0 / 3.0));
  num b = expm1(-1 / (4 * a * a));
  num c = 2 * a * a * b + 0.5;
  return -0.5 * na43 * K *
         (1 - alpha -
          beta * 8.0 / 3.0 * a * (sqrt(M_PI) * erf(1 / (2 * a)) + 2 * a * (b - c)));
}

template <typename num> static num beckex(const densvars<num> & d) {
  return becke_alpha(d.a, d.gaa) + becke_alpha(d.b, d.gbb);
}

template <typename num> static num beckexcorr(const densvars<num> & d) {
  return becke_corr(d.a, d.gaa) + becke_corr(d.b, d.gbb);
}

template <typename num> static num beckesrx(const densvars<num> & d) {
  parameter mu = d.get_param(XC_RANGESEP_MU);
  return becke_sr(mu, d.a, d.gaa) + becke_sr(mu, d.b, d.gbb);
}

template <typename num> static num beckecamx(const densvars<num> & d) {
  parameter mu = d.get_param(XC_RANGESEP_MU);
  parameter alpha = d.get_param(XC_CAM_ALPHA);
  parameter beta = d.get_param(XC_CAM_BETA);
  return becke_cam(alpha, beta, mu, d.a, d.gaa) +
         becke_cam(alpha, beta, mu, d.b, d.gbb);
}

FUNCTIONAL(XC_BECKEX) = {
    "Becke 88 exchange",
    "Becke 88 exchange including Slater part\n"
    "A.D. Becke, Density-functional exchange-energy approximation\n"
    "with correct asymptotic behaviour, Phys. Rev. A38 (1988) 3098-3100.\n"
    "Implemented by Ulf Ekstrom\n"
    "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(beckex) XC_A_B_GAA_GAB_GBB,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-11,
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06, 0.82E+06},
    {-0.277987329958E+03,

     -0.385951846654E+01, -0.381309494319E+01, -0.172434478018E-04,
     0.000000000000E+00,  -0.173712338362E-04,

     -0.441426807406E-01, 0.000000000000E+00,  0.201415922856E-06,
     0.000000000000E+00,  0.000000000000E+00,  -0.447245742260E-01,
     0.000000000000E+00,  0.000000000000E+00,  0.195961359539E-06,
     0.700742719647E-11,  0.000000000000E+00,  0.000000000000E+00,
     0.000000000000E+00,  0.000000000000E+00,  0.718678968862E-11}};

FUNCTIONAL(XC_BECKECORRX) = {
    "Becke 88 exchange correction",
    "Becke 88 exchange not including Slater part\n"
    "A.D. Becke, Density-functional exchange-energy approximation\n"
    "with correct asymptotic behaviour, Phys. Rev. A38 (1988) 3098-3100.\n"
    "Implemented by Ulf Ekstrom\n"
    "Test case from http://www.cse.scitech.ac.uk/ccg/dft/data_pt_x_lda.html\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(beckexcorr) XC_A_B_GAA_GAB_GBB,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-11,
    {0.39e+02, 0.38e+02, 0.81e+06, 0.82e+06, 0.82e+06},
    {
        //     radovan: reference data obtained from *.c implementation in DIRAC
        -3.603918211981e+01, // 00000

        3.479609002901e-01,  // 10000
        3.581112448092e-01,  // 01000
        -1.724344780183e-05, // 00100
        0.000000000000e+00,  // 00001
        -1.737123383621e-05, // 00010

        -8.181318630937e-03,    // 20000
        0.000000000000e+00,     // 11000
        2.014159228564e-07,     // 10100
        0.000000000000e+00,     // 10010
        0.000000000000e+00,     // 10001
        -8.135046261131e-03,    // 02000
        0.000000000000e+00,     // 01100
        0.000000000000e+00,     // 0100
        1.959613595393e-07,     // 01010
        7.0074271964711398e-12, // radovan: i got this using xcfun
        0.0000000000000000,     0.0000000000000000,
        0.0000000000000000,     0.0000000000000000,
        7.1867896886212297e-12, // radovan: i got this using xcfun
    }};

FUNCTIONAL(XC_BECKESRX) = {
    "Short range Becke 88 exchange",
    "Short range Becke 88 exchange, Implemented by Ulf Ekstrom\n"
    "Uses XC_RANGESEP_MU\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(beckesrx)};

FUNCTIONAL(XC_BECKECAMX) = {"CAM Becke 88 exchange",
                            "CAM Becke 88 exchange, Implemented by Elisa Rebolini\n"
                            "Uses XC_RANGESEP_MU\n",
                            XC_DENSITY | XC_GRADIENT,
                            ENERGY_FUNCTION(beckecamx)};
