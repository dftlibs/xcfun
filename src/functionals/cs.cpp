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

template <typename num> static num csc(const densvars<num> & d) {
  parameter a = 1.0;
  parameter b = 1.0;
  parameter c = 1.0;
  parameter dpar = 1.0;
  num gamma = 2 * (1 - (d.a * d.a + d.b * d.b) / (d.n * d.n));
  num curv = d.a * d.taua + d.b * d.taub - (1.0 / 8.0) * d.gnn - (d.jpaa + d.jpbb);
  return -a * gamma *
         (d.n + 2 * b * pow(d.n, -5.0 / 3.0) * curv * exp(-c * d.n_m13)) /
         (1 + dpar * d.n_m13);
}

FUNCTIONAL(XC_CSC) = {"Colle-Salvetti correlation functional",
                      "C. Lee, W. Yang, and R.G. Parr, Development of the \n"
                      "Colle-Salvetti correlation-energy formula into a functional\n"
                      "of the electron density, Phys. Rev. B37 (1988) 785-789\n"
                      "Implemented by Ulf Ekstrom\n",
                      XC_DENSITY | XC_GRADIENT | XC_KINETIC | XC_LAPLACIAN | XC_JP,
                      ENERGY_FUNCTION(csc)};
