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
#include "m0xy_fun.hpp"
#include "pbex.hpp"

// M06 exchange functional. to be used with HF exchange factor of .27

template <typename num> static num m06x(const densvars<num> & d) {
  using m0xy_metagga_xc_internal::alpha_x;
  using m0xy_metagga_xc_internal::fw;
  using m0xy_metagga_xc_internal::h;
  using m0xy_metagga_xc_internal::zet;
  using pw91_like_x_internal::chi2;

  const parameter param_a[12] = {5.877943e-01,
                                 -1.371776e-01,
                                 2.682367e-01,
                                 -2.515898e+00,
                                 -2.978892e+00,
                                 8.710679e+00,
                                 1.688195e+01,
                                 -4.489724e+00,
                                 -3.299983e+01,
                                 -1.449050e+01,
                                 2.043747e+01,
                                 1.256504e+01};
  const parameter param_d[6] = {1.422057e-01,
                                7.370319e-04,
                                -1.601373e-02,
                                0.000000e+00,
                                0.000000e+00,
                                0.000000e+00};

  num chia2 = chi2(d.a, d.gaa);
  num chib2 = chi2(d.b, d.gbb);

  return ((pbex::energy_pbe_ab(pbex::R_pbe, d.a, d.gaa) * fw(param_a, d.a, d.taua) +
           lsda_x(d.a) * h(param_d, alpha_x, chia2, zet(d.a, d.taua))) +
          (pbex::energy_pbe_ab(pbex::R_pbe, d.b, d.gbb) * fw(param_a, d.b, d.taub) +
           lsda_x(d.b) * h(param_d, alpha_x, chib2, zet(d.b, d.taub))));
}

FUNCTIONAL(XC_M06X) = {
    "M06 exchange",
    "M06 Meta-Hybrid Exchange Functional\n"
    "Y Zhao and D. G. Truhlar, Theor. Chem. Account 120, 215 (2008)\n"
    "Implemented by Andre Gomes\n"
    "Reference gradient from ADF\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC,
    ENERGY_FUNCTION(m06x)
#ifdef XCFUN_REF_PBEX_MU
        XC_A_B_GAA_GAB_GBB_TAUA_TAUB,
    XC_PARTIAL_DERIVATIVES,
    1,
    1e-7,
    {0.153652558932587,
     0.153652558932587,
     8.390981882024769E-002,
     8.390981882024769E-002,
     8.390981882024769E-002,
     6.826262722466429E-002,
     6.826262722466429E-002},
    {-1.0989923683183833e-01, // Energy from xcfun
     -0.519568215542419,
     -0.519568215542419,
     -1.757089749188437E-002,
     0.000000000000000E+000,
     -1.757089749188437E-002,
     9.227728385689479E-002,
     9.227728385689479E-002}
#endif
};

/*
const double d[] =
  {1., .8, 1., 1., 1., .33, .21};
const double ref[] =
  { -0.93035619, -1.57676773, -1.59326350, -0.00081086,   0.00000000, -0.00129163,
1.62369085,  2.06708653 };
*/
