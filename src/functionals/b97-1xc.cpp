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

template <typename num> static num b97_1x_en(const densvars<num> & d) {
  return b97x::energy_b97x_ab(b97x::Gamma, b97x::c_b97_1, d.a_43, d.gaa) +
         b97x::energy_b97x_ab(b97x::Gamma, b97x::c_b97_1, d.b_43, d.gbb);
}

template <typename num> static num b97_1c_en(const densvars<num> & d) {
  num e_LSDA_a, e_LSDA_b, tmp;
  tmp = b97c::energy_b97c_par(
            b97c::Gamma_par, b97c::c_b97_1[1], d.a, d.a_43, d.gaa, e_LSDA_a) +
        b97c::energy_b97c_par(
            b97c::Gamma_par, b97c::c_b97_1[1], d.b, d.b_43, d.gbb, e_LSDA_b);

  return tmp + b97c::energy_b97c_antipar(
                   b97c::Gamma_antipar, b97c::c_b97_1[0], d, e_LSDA_a, e_LSDA_b);
}

FUNCTIONAL(XC_B97_1X) = {"B97-1 exchange",
                         "Hybrid exchange-correlation functional determined from\n"
                         "thermochemical data NO ab initio potentials\n"
                         "HCTH; J.Chem.Phys.; 109, 6264,  (1998)\n"
                         "Implemented by Alex Borgoo ",
                         XC_DENSITY | XC_GRADIENT,
                         ENERGY_FUNCTION(b97_1x_en)};

FUNCTIONAL(XC_B97_1C) = {"B97-1 correlation",
                         "Hybrid exchange-correlation functional determined from\n"
                         "thermochemical data NO ab initio potentials\n"
                         "HCTH; J.Chem.Phys.; 109, 6264, (1998)\n"
                         "Implemented by Alex Borgoo ",
                         XC_DENSITY | XC_GRADIENT,
                         ENERGY_FUNCTION(b97_1c_en)};
