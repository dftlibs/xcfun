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

template <typename num> static num pbeloc_eps(const densvars<num> & d) {
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
  return (eps + H);
}

template <typename num> static num pbeloc_eps_pola(const num & a, const num & gaa) {
  using xc_constants::param_gamma;
  const parameter beta0 = 0.0375;
  const parameter aa = 0.08;
  num u = pow(2.0, -1.0 / 3.0); // phi for fully polarized systems
  num u3 = pow3(u);
  // d2 is t^2
  num d2 = pow(1.0 / 12 * pow(3, 5.0 / 6.0) / pow(M_PI, -1.0 / 6), 2.0) * gaa /
           (u * u * pow(a, 7.0 / 3.0));
  num rs = pow(3 / (4 * M_PI), 1.0 / 3.0) * pow(a, -1.0 / 3.0);
  num ff = 1 - exp(-rs * rs);
  num beta = beta0 + aa * d2 * ff;
  num bg = beta / param_gamma;
  num eps = pw92eps::pw92eps_polarized(a);
  num A = bg / expm1(-eps / (param_gamma * u3));
  num d2A = d2 * A;
  num H = param_gamma * u3 * log(1 + bg * d2 * (1 + d2A) / (1 + d2A * (1 + d2A)));
  return (eps + H);
}

template <typename num> static num C(const densvars<num> & d) {
  num gzeta2 = (pow2(d.n) * d.gss - 2 * d.n * d.s * d.gns + pow2(d.s) * d.gnn) /
               pow(d.n, 4);                                      // (grad zeta)^2
  num xi2 = gzeta2 / (4 * pow(3 * pow2(M_PI) * d.n, 2.0 / 3.0)); // xi^2
  num C0 = 0.35 + 0.87 * pow2(d.zeta) + 0.50 * pow(d.zeta, 4) +
           2.26 * pow(d.zeta, 6); // C(zeta,0)
  return C0 * pow(1 + 0.5 * xi2 * (ufunc(d.zeta, -4.0 / 3.0)), -4);
}

template <typename num>
static num epsc_summax(const densvars<num> & d)
// sum of mmax between fully and partially polarized PBEloc energies per particle
{
  num epsc_pbeloc = pbeloc_eps(d);
  num epsc_pbeloc_a = pbeloc_eps_pola(d.a, d.gaa);
  num epsc_pbeloc_b = pbeloc_eps_pola(d.b, d.gbb);
  return (d.a * max(epsc_pbeloc, epsc_pbeloc_a) +
          d.b * max(epsc_pbeloc, epsc_pbeloc_b)) /
         d.n;
}

template <typename num> static num epsc_revpkzb(const densvars<num> & d) {
  num tauwtau2 = pow2(d.gnn / (8.0 * d.n * d.tau));
  num epsc_sum = epsc_summax(d);
  num epsc_pbeloc = pbeloc_eps(d);
  num CC = C(d);
  return epsc_pbeloc * (1 + CC * tauwtau2) - (1 + CC) * tauwtau2 * epsc_sum;
}

template <typename num> static num energy(const densvars<num> & d) {
  num eps_pkzb = epsc_revpkzb(d);
  num tauwtau3 = pow3(d.gnn / (8.0 * d.n * d.tau));
  const parameter dd = 4.5;
  return d.n * eps_pkzb * (1 + dd * eps_pkzb * tauwtau3);
}

FUNCTIONAL(XC_TPSSLOCC) = {"TPSSloc correlation functional",
                           "TPSSloc correlation functional.\n"
                           "L.A. Constantin, E.Fabiano, F. Della Sala,\n"
                           "      Phys. Rev. B 86, 035130 (2012).\n"
                           "Implemented by Eduardo Fabiano\n",
                           XC_DENSITY | XC_GRADIENT | XC_KINETIC,
                           ENERGY_FUNCTION(energy)};
