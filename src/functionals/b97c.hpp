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

#include "b97xc.hpp"
#include "functional.hpp"
#include "pw92eps.hpp"

namespace b97c {

// parameters for enhancement factor
// c_ab0, c_ab1, c_ab2, c_ss0, c_ss1, c_ss2
const parameter c_b97[2][3] = {{0.9454, 0.7471, -4.5961}, {0.1737, 2.3487, -2.4868}};

const parameter c_b97_1[2][3] = {{0.955689, 0.788552, -5.47869},
                                 {0.0820011, 2.71681, -2.87103}};

const parameter c_b97_2[2][3] = {{0.999849, 1.40626, -7.44060},
                                 {0.585808, -0.691682, 0.394796}};

const parameter Gamma_par = 0.2;
const parameter Gamma_antipar = 0.006;

template <typename num>
static num energy_b97c_par(const parameter & Gamma,
                           const parameter c_params[],
                           const num & a,
                           const num & a_43,
                           const num & gaa,
                           num & e_LSDA) {
  e_LSDA = pw92eps::pw92eps_polarized(a) * a;
  num s2_ab2 = b97xc::spin_dens_gradient_ab2(gaa, a_43);
  return e_LSDA * b97xc::enhancement(Gamma, c_params, s2_ab2);
}

template <typename num>
static num energy_b97c_antipar(const parameter & Gamma,
                               const parameter c_params[],
                               const densvars<num> & d,
                               const num & e_LSDA_a,
                               const num & e_LSDA_b) {
  num e_LSDA = pw92eps::pw92eps(d) * d.n - e_LSDA_a - e_LSDA_b;
  num s2_ab = 0.5 * (b97xc::spin_dens_gradient_ab2(d.gaa, d.a_43) +
                     b97xc::spin_dens_gradient_ab2(d.gbb, d.b_43));

  return e_LSDA * b97xc::enhancement(Gamma, c_params, s2_ab);
}
} // namespace b97c
