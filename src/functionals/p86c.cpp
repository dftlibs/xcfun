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
#include "pz81c.hpp"

template <typename num> static num Cg(const num & r) {
  parameter Cx = 0.001667;
  parameter Bg = 0.000007389;
  return Cx + (0.002568 + r * (0.023266 + Bg * r)) /
                  (1 + r * (8.723 + r * (0.472 + 10000 * Bg * r)));
}

template <typename num> static num Pg(const densvars<num> & d) {
  parameter Fg = 0.11;
  parameter Cinf = 0.004235;
  parameter fudge = 1e-12; // Avoid instability at d.gnn = 0
#ifndef INEXACT_PI
  parameter pi_expr = pow(9 * M_PI, 1.0 / 6.0);
#else
  parameter pi_expr = 1.745;
#endif
  return pi_expr * Fg * Cinf * sqrt(fudge + d.gnn) /
         (Cg(d.r_s) * pow(d.n, 7.0 / 6.0));
}

template <typename num> static num dz(const densvars<num> & d) {
  return cbrt(2.0) * sqrt(pow(d.a, 5.0 / 3.0) + pow(d.b, 5.0 / 3.0)) *
         pow(d.n, -5.0 / 6.0);
}

template <typename num> static num p86c(const densvars<num> & d) {
  return d.n * pz81eps::pz81eps(d) +
         exp(-Pg(d)) * Cg(d.r_s) * d.gnn / (pow(d.n, 4.0 / 3.0) * dz(d));
}

template <typename num> static num p86c_corr(const densvars<num> & d) {
  return exp(-Pg(d)) * Cg(d.r_s) * d.gnn / (pow(d.n, 4.0 / 3.0) * dz(d));
}

FUNCTIONAL(XC_P86C) = {
    "P86C GGA correlation",
    "J.P. Density-functional approximation for the correlation energy\n"
    "of the inhomogeneous electron , Phys. Rev. B, 33(12):8822gasPerdew,\n"
    "Implemented by Ulf Ekstrom.\n"
    "Reference data from ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_c_p86.html",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(p86c) XC_A_B_GAA_GAB_GBB,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-10,
#ifdef HIGH_DENSITY
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06, 0.82E+06},
#ifdef INEXACT_PI
    {
        -0.356963343227E+01, -0.433530899660E-01, -0.447602011737E-01,
        -0.122488989075E-06, -0.244977978150E-06, -0.122488989075E-06,
        -0.169450157064E-02, -0.308422706942E-02, 0.227151852191E-07,
        0.454303704383E-07,  0.227151852191E-07,  -0.165896811718E-02,
        0.226922323054E-07,  0.453844646108E-07,  0.226922323054E-07,
        -0.207808284633E-12, -0.415616569266E-12, -0.207808284633E-12,
        -0.831233138532E-12, -0.415616569266E-12, -0.207808284633E-12,
    }
#else
    {
        -3.5721593740706306e+00, -4.3353602806095654e-02, -4.4761187343925160e-02,
        -1.2281210551047686e-07, -2.4562421102095372e-07, -1.2281210551047686e-07,
        -1.6930136257204879e-03, -3.0832001805399236e-03, +2.2702465560744966e-08,
        +4.5404931121489932e-08, +2.2702465560744966e-08, -1.6574558665771906e-03,
        +2.2679452099010038e-08, +4.5358904198020076e-08, +2.2679452099010038e-08,
        -2.0762304168807950e-13, -4.1524608337615900e-13, -2.0762304168807950e-13,
        -8.3049216675231800e-13, -4.1524608337615900e-13, -2.0762304168807950e-13,
    }
#endif
#else
    {0.48E-01, 0.25E-01, 0.46E-02, 0.44E-02, 0.41E-02},
#ifdef INEXACT_PI
    {-0.227465804486E-02, -0.607101652780E-01, -0.809996785446E-01,
     0.542965890398E-01,  0.108593178080E+00,  0.542965890398E-01,
     0.329942643730E+00,  -0.517763280120E+00, -0.232199446020E+00,
     -0.464398892040E+00, -0.232199446020E+00, 0.564642641219E+00,
     0.164136266720E-01,  0.328272533439E-01,  0.164136266720E-01,
     -0.145741442973E+01, -0.291482885947E+01, -0.145741442973E+01,
     -0.582965771894E+01, -0.291482885947E+01, -0.145741442973E+01}
#else // self-computed values
    {
        -2.2748317487871362e-03, -6.0705512544759296e-02, -8.0995821166401388e-02,
        +5.4284455865926534e-02, +1.0856891173185307e-01, +5.4284455865926534e-02,
        +3.2978359263409207e-01, -5.1791169779733182e-01, -2.3191926560979093e-01,
        -4.6383853121958185e-01, -2.3191926560979093e-01, +5.6456418726387536e-01,
        +1.6638251739185206e-02, +3.3276503478370412e-02, +1.6638251739185206e-02,
        -1.4574898694069474e+00, -2.9149797388138947e+00, -1.4574898694069474e+00,
        -5.8299594776277894e+00, -2.9149797388138947e+00, -1.4574898694069474e+00,
    }
#endif
#endif
};

FUNCTIONAL(XC_P86CORRC) = {
    "P86C GGA correlation",
    "J.P. Density-functional approximation for the correlation energy\n"
    "of the inhomogeneous electron , Phys. Rev. B, 33(12):8822gasPerdew,\n"
    "Implemented by Ulf Ekstrom.\n"
    "Reference data from ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_c_p86.html",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(p86c_corr)};
