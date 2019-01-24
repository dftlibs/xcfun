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

#include "b97c.hpp"
#include "b97x.hpp"
#include "constants.hpp"
#include "functional.hpp"

template <typename num> static num b97x_en(const densvars<num> & d) {
  return b97x::energy_b97x_ab(b97x::Gamma, b97x::c_b97, d.a_43, d.gaa) +
         b97x::energy_b97x_ab(b97x::Gamma, b97x::c_b97, d.b_43, d.gbb);
}

template <typename num> static num b97c_en(const densvars<num> & d) {
  num e_LSDA_a, e_LSDA_b, tmp;
  tmp = b97c::energy_b97c_par(
            b97c::Gamma_par, b97c::c_b97[1], d.a, d.a_43, d.gaa, e_LSDA_a) +
        b97c::energy_b97c_par(
            b97c::Gamma_par, b97c::c_b97[1], d.b, d.b_43, d.gbb, e_LSDA_b);

  return tmp + b97c::energy_b97c_antipar(
                   b97c::Gamma_antipar, b97c::c_b97[0], d, e_LSDA_a, e_LSDA_b);
}

FUNCTIONAL(XC_B97X) = {
    "B97 exchange",
    "Density-functional thermochemistry. V.\n"
    "Systematic optimization of exchange-correlation functionals\n"
    "A.D.Becke; J.Chem.Phys.; 107, 8554, (1997)\n"
    "Implemented by Sarah Reimann ",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(b97x_en)};

FUNCTIONAL(XC_B97C) = {
    "B97 correlation",
    "Density-functional thermochemistry\n"
    "Systematic optimization of exchange-correlation functionals\n"
    "A.D.Becke; J.Chem.Phys.; 107, 8554, (1997)\n"
    "Implemented by Sarah Reimann ",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(b97c_en)};
