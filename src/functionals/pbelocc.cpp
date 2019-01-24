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
#include "pw92eps.hpp"

// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template <typename num> static num phi(const densvars<num> & d) {
  return pow(2.0, -1.0 / 3.0) * d.n_m13 * d.n_m13 * (sqrt(d.a_43) + sqrt(d.b_43));
}

template <typename num> static num energy(const densvars<num> & d) {
  using xc_constants::param_gamma;
  const parameter beta0 = 0.0375;
  const parameter aa = 0.08;
  num u = phi(d);
  num u3 = pow3(u);
  // d2 is t^2
  num d2 = pow(1.0 / 12 * pow(3, 5.0 / 6.0) / pow(M_PI, -1.0 / 6), 2.0) * d.gnn /
           (u * u * pow(d.n, 7.0 / 3.0));
  num ff = 1 - exp(-d.r_s * d.r_s);
  num beta = beta0 + aa * d2 * ff;
  num bg = beta / param_gamma;
  num eps = pw92eps::pw92eps(d);
  num A = bg / expm1(-eps / (param_gamma * u3));
  num d2A = d2 * A;
  num H = param_gamma * u3 * log(1 + bg * d2 * (1 + d2A) / (1 + d2A * (1 + d2A)));
  return d.n * (eps + H);
}

FUNCTIONAL(XC_PBELOCC) = {
    "PBEloc correlation functional.",
    "PBEloc correlation functional\n"
    "L.A. Constantin, E. Fabiano, S. Laricchia, F. Della Sala,\n"
    "      Phys. Rev. B 86, 035130 (2012).\n"
    "Implemented by Eduardo Fabiano\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(energy)};
