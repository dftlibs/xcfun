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

#include "constants.hpp"
#include "functional.hpp"

template <typename num>
static num energy_blocx(const num & d_n, const num & d_gnn, const num & d_tau) {
  const parameter kappa = 0.804;
  const parameter mu = 0.21951;
  const parameter b = 0.40;
  const parameter e = 1.537;
  const parameter c = 1.59096;
  num p0 = 1.0 / (4 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
  num p = d_gnn * p0; // s^2
  num tauw = d_gnn / (8.0 * d_n);
  num z = tauw / d_tau;
  num z2 = pow2(z);
  num tau_unif = 0.3 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 5.0 / 3.0);
  num alpha = (d_tau - tauw) / tau_unif; //  5*p*z/(3.0 - 3.0*z);
  num q_b =
      (9.0 / 20.0) * (alpha - 1) / sqrt(1 + b * alpha * (alpha - 1)) + 2 * p / 3.0;
  num ff = 4 - 3.3 * z;
  num zf = exp(log(z) * ff); // z^f
  num tmp1 = p * (10.0 / 81.0 + c * zf / pow2(1 + z2));
  num tmp2 = 146.0 * pow2(q_b) / 2025.0;
  num tmp3 = -73.0 / 405.0 * q_b * d_gnn *
             sqrt((0.5 * 0.6 * 0.6) * pow(8 * d_n * d_tau, -2) + 0.5 * p0 * p0);
  num tmp4 = pow2(p * 10.0 / 81.0) / kappa;
  num tmp5 = 2 * sqrt(e) * 0.6 * 0.6 * z2 * 10.0 / 81.0 + e * mu * pow3(p);
  num tmp6 = tmp1 + tmp2 + tmp3 + tmp4 + tmp5;
  num x = tmp6 / pow2(1 + sqrt(e) * p);
  num Fx = 1 + kappa - kappa / (1 + x / kappa);
  num lda = (-0.75 * pow(3 / PI, 1.0 / 3.0)) * pow(d_n, 4.0 / 3.0); // lda exchange
  return lda * Fx;
}

template <typename num> static num energy(const densvars<num> & d) {
  num enea = energy_blocx(2 * d.a, 4 * d.gaa, 2 * d.taua);
  num eneb = energy_blocx(2 * d.b, 4 * d.gbb, 2 * d.taub);
  return (enea + eneb) / 2;
}

FUNCTIONAL(XC_BLOCX) = {"BLOC exchange functional",
                        "BLOC exchange functional.\n"
                        "L.A. Constantin, E.Fabiano, F. Della Sala,\n"
                        "      J. Chem. Theory Comput. 9, 2256 (2013).\n"
                        "Implemented by Eduardo Fabiano\n",
                        XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                        ENERGY_FUNCTION(energy)};
