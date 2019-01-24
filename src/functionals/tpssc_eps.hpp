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

#pragma once

#include "constants.hpp"
#include "functional.hpp"
#include "pbec_eps.hpp"

namespace tpssc_eps {
template <typename num> static num C(const densvars<num> & d) {
  num gzeta2 = (pow2(d.n) * d.gss - 2 * d.n * d.s * d.gns + pow2(d.s) * d.gnn) /
               pow(d.n, 4);                                      // (grad zeta)^2
  num xi2 = gzeta2 / (4 * pow(3 * pow2(M_PI) * d.n, 2.0 / 3.0)); // xi^2
  num C0 = 0.53 + 0.87 * pow2(d.zeta) + 0.50 * pow(d.zeta, 4) +
           2.26 * pow(d.zeta, 6); // C(zeta,0)
  return C0 * pow(1 + 0.5 * xi2 * (ufunc(d.zeta, -4.0 / 3.0)), -4);
  //    = return C0/pow( 1 + 0.5*xi2*(pow(1+d.zeta,-4.0/3.0) +
  //    pow(1-d.zeta,-4.0/3.0)),4);
}

template <typename num>
static num epsc_summax(const densvars<num> & d) /* sum_sigma n_sigma*epsc_max/n is
                                                   the sum in eq [12] of the
                                                   reference */
{
  num epsc_pbe = pbec_eps::pbec_eps(d);
  num epsc_pbe_a = pbec_eps::pbec_eps_polarized(d.a, d.gaa);
  num epsc_pbe_b = pbec_eps::pbec_eps_polarized(d.b, d.gbb);
  return (d.a * max(epsc_pbe, epsc_pbe_a) + d.b * max(epsc_pbe, epsc_pbe_b)) / d.n;
}

template <typename num>
static num epsc_revpkzb(const densvars<num> & d) /* eps_c^{rev PKZB} it is the eq
                                                    [12] from the reference */
{
  num tauwtau2 = pow2(d.gnn / (8.0 * d.n * d.tau));
  num epsc_sum = epsc_summax(d);
  num epsc_pbe = pbec_eps::pbec_eps(d);
  num C_zeta_xi = C(d);
  return epsc_pbe * (1 + C_zeta_xi * tauwtau2) -
         (1 + C_zeta_xi) * tauwtau2 * epsc_sum;
}

template <typename num> static num tpssc_eps(const densvars<num> & d) {
  num eps_pkzb = epsc_revpkzb(d);
  num tauwtau3 = pow3(d.gnn / (8.0 * d.n * d.tau));
  parameter dd = 2.8;
  return eps_pkzb * (1 + dd * eps_pkzb * tauwtau3);
}
} // namespace tpssc_eps
