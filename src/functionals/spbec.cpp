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
#include "vwn.hpp"

/* "Simplified" PBEc, spbec. Implemented by ulfek */

/* These are the constants in the paper, not very high precision.. */
static const parameter gamm = 0.066725;
static const parameter beta = 0.031091;
static const parameter beta_gamma = beta / gamm;

template <typename num> static num G(const num & eps, const num & phi3) {
  return beta_gamma / expm1(-eps / (beta_gamma * phi3));
}

template <typename num>
static num H_spbe(const num & t2, const num & eps, const num & phi3) {
  return gamm * phi3 * log(1 + beta_gamma * t2 / (1 + t2 * G(eps, phi3)));
}

// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template <typename num> static num phi(const densvars<num> & d) {
  return pow(2.0, -1.0 / 3.0) * d.n_m13 * d.n_m13 * (sqrt(d.a_43) + sqrt(d.b_43));
}

template <typename num> static num spbec(const densvars<num> & d) {
  num eps = vwn::vwn5_eps(d);
  num p = phi(d);
  num t2 = (cbrt(M_PI / 3) / 16) * d.gnn * d.n_m13 / pow2(p * d.n);
  return d.n * (eps + H_spbe(t2, eps, pow3(p)));
}

FUNCTIONAL(XC_SPBEC) = {
    "sPBE correlation functional",
    "Simplified PBE correlation functional for use with the SSB functionals.\n"
    "Swart, M. and Sola, M. and Bickelhaupt M.; JCP 131 094103 (2009)\n"
    "Implemented by Ulf Ekstrom\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(spbec)};
