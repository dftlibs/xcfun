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
#include "pw92eps.hpp"

/*
  This code was adapted by Ulf Ekstrom from a Fortran program provided
  by Paola Gori-Giorgi.
 */

template <typename num> static num Qrpa(const num & x) {
  const parameter Acoul = 2.0 * (log(2.0) - 1.0) / (M_PI * M_PI);
  const parameter a2 = 5.84605;
  const parameter c2 = 3.91744;
  const parameter d2 = 3.44851;
  const parameter b2 =
      d2 - 3.0 / (2 * M_PI * Acoul) * pow(4.0 / (9.0 * M_PI), 1.0 / 3.0);
  return Acoul * log((1 + x * (a2 + x * (b2 + c2 * x))) / (1 + x * (a2 + d2 * x)));
}

template <typename num> static num dpol(const num & rs) {
  const parameter cf = pow(9.0 * M_PI / 4.0, 1.0 / 3.0);
  const parameter p2p = 0.04;
  const parameter p3p = 0.4319;
  num rs2 = rs * rs;
  return pow(2.0, 5.0 / 3.0) / 5.0 * pow(cf, 2) / rs2 *
         (1.0 + (p3p - 0.454555) * rs) / (1.0 + p3p * rs + p2p * rs2);
}
/*
 on-top pair-distribution function
 Gori-Giorgi and Perdew, PRB 64, 155102 (2001)
 x -> rs
*/
template <typename num> static num g0f(const num & x) {
  const parameter C0f = 0.0819306;
  const parameter D0f = 0.752411;
  const parameter E0f = -0.0127713;
  const parameter F0f = 0.00185898;
  return (1 + x * (D0f - 0.7317 + x * (C0f + x * (E0f + F0f * x)))) * exp(-D0f * x) /
         2.0;
}

template <typename num>
static num ecorrlr(const densvars<num> & d, parameter mu, const num & ec) {
  const parameter alpha = pow(4.0 / 9.0 / M_PI, 1.0 / 3.0);
  const parameter cf = 1 / alpha; // TODO: The normal CF?
  num phi = (pow(1.0 + d.zeta, 2.0 / 3.0) + pow(1.0 - d.zeta, 2.0 / 3.0)) / 2.0;
  // cc parameters from the fit
  const parameter adib = 0.784949;
  const parameter q1a = -0.388;
  const parameter q2a = 0.676;
  const parameter q3a = 0.547;
  const parameter t1a = -4.95;
  const parameter t2a = 1.0;
  const parameter t3a = 0.31;

  num b0 = adib * d.r_s;
  num rs2 = d.r_s * d.r_s;
  num rs3 = rs2 * d.r_s;

  num d2anti = (q1a * d.r_s + q2a * rs2) * exp(-q3a * d.r_s) / rs2;
  num d3anti = (t1a * d.r_s + t2a * rs2) * exp(-t3a * d.r_s) / rs3;

  const num & z = d.zeta;
  num z2 = d.zeta * d.zeta;
  num coe2 = -3.0 / 8.0 / rs3 * (1.0 - z2) * (g0f(d.r_s) - 0.5);

  num coe3 = -(1.0 - z2) * g0f(d.r_s) / (sqrt(2.0 * M_PI) * rs3);

  num coe4 =
      -9.0 / 64.0 / rs3 *
      (pow((1.0 + z) / 2.0, 2) * dpol(d.r_s * pow(2 / (1.0 + z), 1.0 / 3.0)) +
       pow((1.0 - z) / 2.0, 2) * dpol(d.r_s * pow(2.0 / (1.0 - z), 1.0 / 3.0)) +
       (1 - z * z) * d2anti -
       pow(cf, 2) / 10.0 * (pow(1.0 + z, 8.0 / 3.0) + pow(1 - z, 8.0 / 3.0)) / rs2);
  num coe5 =
      -9.0 / 40.0 / (sqrt(2.0 * M_PI) * rs3) *
      (pow((1.0 + z) / 2.0, 2) * dpol(d.r_s * pow(2.0 / (1.0 + z), 1.0 / 3.0)) +
       pow((1.0 - z) / 2.0, 2) * dpol(d.r_s * pow(2.0 / (1.0 - z), 1.0 / 3.0)) +
       (1.0 - z2) * d3anti);

  num b06 = pow(b0, 6);
  num b08 = pow(b0, 8);
  num a1 = 4 * b06 * coe3 + b08 * coe5;
  num a2 = 4 * b06 * coe2 + b08 * coe4 + 6 * pow(b0, 4) * ec;
  num a3 = b08 * coe3;
  num a4 = b06 * (pow(b0, 2) * coe2 + 4 * ec);

  return (pow(phi, 3) * Qrpa(mu * sqrt(d.r_s) / phi) + a1 * pow(mu, 3) +
          a2 * pow(mu, 4) + a3 * pow(mu, 5) + a4 * pow(mu, 6) +
          pow(b0 * mu, 8) * ec) /
         pow(1 + pow(b0 * mu, 2), 4);
}
template <typename num> static num ldaerfc(const densvars<num> & d) {
  double mu = d.get_param(XC_RANGESEP_MU);
  num eps = pw92eps::pw92eps(d);
  return d.n * (eps - ecorrlr(d, mu, eps));
}

FUNCTIONAL(XC_LDAERFC) = {
    "Short-range spin-dependent LDA correlation functional",
    "Short-range spin-dependent LDA correlation functional from\n"
    "Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)"
    "Adapted from Gori-Giorgi and MOLPRO by Ulf Ekstrom\n"
    "Test case from Gori-Giorgi (personal communication),\n"
    "up to 10^-7, then xcfun decimals due to more accurate pw92c.\n"
    "Range separation parameter is XC_RANGESEP_MU\n",
    XC_DENSITY,
    ENERGY_FUNCTION(ldaerfc) XC_A_B,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-7,
    {1.1, 1.0},
    {
        -1.4579390272267870e-01, // Generated with XC_RANGESEP_MU=0.4
        -7.7624817385549980e-02,
        -8.2132052511772885e-02,
        +1.5795011054215363e-02,
        -2.7440928179985190e-02,
        +1.9539616096309973e-02,
    }};
/* Original numbers from Paola, with inaccurate pw92c
const double ref[] =
  {
    -0.1457945300494694,
    -0.07762517247521351,
    -0.08213242046922536,
    0.015795038615569225,
    -0.02744102459091926,
    0.019539653410626807
  };
*/
