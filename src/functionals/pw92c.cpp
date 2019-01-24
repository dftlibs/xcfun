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

template <typename num> static num pw92c(const densvars<num> & d) {
  return pw92eps::pw92eps(d) * d.n;
}

FUNCTIONAL(XC_PW92C) = {"PW92 LDA correlation",
                        "Accurate and simple analytic representation of the "
                        "electron-gas correlation energy\n"
                        "J.P.Perdew, Y. Wang; Phys. Rev. B; 45, 13244, (1992)\n"
                        "Implemented by Ulf Ekstrom. Some parameters have higher\n"
                        "accuracy than given in the paper.\n",
                        XC_DENSITY,
                        ENERGY_FUNCTION(pw92c) XC_A_B,
                        XC_PARTIAL_DERIVATIVES,
                        2,
                        1e-11,
                        {0.39E+02, 0.38E+02},
                        {
                            -8.4713855882783946e+00, // self computed
                            -1.1861930857502517e-01,
                            -1.2041769989725633e-01,
                            +7.5202855619095870e-04,
                            -1.0249091426230799e-03,
                            +7.9516089195232130e-04,
                        }};

/* Daresbury numbers:
const double ref[] =
  {
    -0.847142825874E+01,
    -0.118619938062E+00,
    -0.120418335387E+00,
    0.752030427237E-03,
    -0.102491320671E-02,
    0.795162900251E-03,
  };
*/
