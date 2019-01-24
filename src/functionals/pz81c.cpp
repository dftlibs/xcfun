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

#include "pz81c.hpp"
#include "functional.hpp"

template <typename num> static num pz81c(const densvars<num> & d) {
  return pz81eps::pz81eps(d) * d.n;
}

FUNCTIONAL(XC_PZ81C) = {"PZ81 LDA correlation",
                        "Implemented by Ulf Ekstrom. Test from "
                        "http://www.cse.scitech.ac.uk/ccg/dft/data_pt_c_pz81.html\n",
                        XC_DENSITY,
                        ENERGY_FUNCTION(pz81c) XC_A_B,
                        XC_PARTIAL_DERIVATIVES,
                        2,
                        1e-11,
#ifdef HIGH_DENSITY
                        {0.39E+02, 0.38E+02},
                        {
                            -0.847966726388E+01,
                            -0.118689256817E+00,
                            -0.121016447321E+00,
                            0.100782705799E-02,
                            -0.129204915717E-02,
                            0.106274642401E-02,
                        }
#else
                        {0.48E-01, 0.25E-01},
                        {-0.358997585489E-02,
                         -0.468661877874E-01,
                         -0.731782746282E-01,
                         0.218577885080E+00,
                         -0.646538277526E+00,
                         0.867717298846E+00}
#endif
};
