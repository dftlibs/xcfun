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

template <typename num>
static num Gc(const num & r,
              const parameter A,
              const parameter a1,
              const parameter b1,
              const parameter b2,
              const parameter b3,
              const parameter b4,
              const parameter p) {
  num sqrtr = sqrt(r);
  return -2 * A * (1 + a1 * r) *
         log(1 + 0.5 / (A * (sqrtr * (b1 + sqrtr * (b2 + sqrtr * b3)) +
                             b4 * pow(r, p + 1))));
}

// (1+(a-b)/(a+b))^p + (1-(a-b)/(a+b))^p =
// ((b^p+a^p)*2^p)/(b+a)^p
template <typename num> static num uf(const densvars<num> & d, const parameter p) {
  return (pow(d.a, p) + pow(d.b, p)) * pow(2 / d.n, p);
}

template <typename num> static num pw91c(const densvars<num> & d) {
  const parameter pa = 1.0;
  const parameter Aa = 0.016887; // ok
  const parameter a1a = 0.11125; // ok
  const parameter b1a = 10.357;  // ok
  const parameter b2a = 3.6231;  // ok
  const parameter b3a = 0.88026; // ok
  const parameter b4a = 0.49671; // ok
  const parameter pe = 1;
  const parameter c0p = 0.031091;   // ok
  const parameter a1p = 0.21370;    // ok
  const parameter b1p = 7.5957;     // ok
  const parameter b2p = 3.5876;     // ok
  const parameter b3p = 1.6382;     // ok
  const parameter b4p = 0.49294;    // ok
  const parameter c0f = 0.015545;   // ok
  const parameter a1f = 0.20548;    // ok
  const parameter b1f = 14.1189;    // ok
  const parameter b2f = 6.1977;     // ok
  const parameter b3f = 3.3662;     // ok
  const parameter b4f = 0.62517;    // ok
  const parameter d2fz0 = 1.709921; // ok
  num fz = (uf(d, 4.0 / 3.0) - 2) / (2 * pow(2, 1.0 / 3.0) - 2);
  num Ac = Gc(d.r_s, Aa, a1a, b1a, b2a, b3a, b4a, pa);
  num EcP = Gc(d.r_s, c0p, a1p, b1p, b2p, b3p, b4p, pe);
  num EcF = Gc(d.r_s, c0f, a1f, b1f, b2f, b3f, b4f, pe);
  num Ec = EcP - Ac * fz * (1 - pow(d.zeta, 4)) / d2fz0 +
           (EcF - EcP) * fz * pow(d.zeta, 4);
  num kF = cbrt(3 * M_PI * M_PI * d.n);
  num ks = sqrt(4) * sqrt(kF / M_PI);
  num gs = 0.5 * uf(d, 2.0 / 3.0);
  num T2 = 0.25 * d.gnn / pow(gs * ks * d.n, 2);
  const parameter alpha = 0.09;
  const parameter Cc0 = 0.004235;
  const parameter Cx = -0.001667;
  const parameter nu = 16 * cbrt(3 * M_PI * M_PI) / M_PI;
  const parameter beta = nu * Cc0;
  num A = 2 * alpha / beta / expm1(-2 * (alpha * Ec) / (pow(gs, 3) * beta * beta));
  num Cc = 1.0 / 1000 *
               ((2.568 + d.r_s * (23.266 + 0.007389 * d.r_s)) /
                (1 + d.r_s * (8.723 + d.r_s * (0.472 + d.r_s * 0.073890)))) -
           Cx;
  num H0 = 0.5 * pow(gs, 3) * beta * beta / alpha *
           log((1 + 2 * (alpha * (T2 + A * T2 * T2)) /
                        (beta * (1 + A * T2 * (1 + A * T2)))));
  num H1 = nu * (Cc - Cc0 - 3.0 / 7.0 * Cx) * pow(gs, 3) * T2 *
           exp(-100 * pow(gs, 4) * pow(ks, 2) * T2 / pow(kF, 2));
  return d.n * (Ec + H0 + H1);
}

#define TEST2
FUNCTIONAL(XC_PW91C) = {
    "PW91 Correlation",
    "J.P. Perdew, J.A. Chevary, S.H. Vosko, K.A. Jackson, M.R. Pederson, D.J. Singh "
    "and"
    "C. Fiolhais, 'Atoms, molecules, solids and surfaces: Applications of the "
    "generalized"
    "gradient approximation for exchange and correlation', Phys. Rev. B, 46(11):6671"
    "6687, 1992"
    "Implemented by Ulf Ekstrom. Test from "
    "ftp://ftp.dl.ac.uk/qcg/dft_library/data_pt_c_pw91.html\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(pw91c) XC_A_B_GAA_GAB_GBB,
    XC_PARTIAL_DERIVATIVES,
    2,     // test order
    1e-11, // test threshold
#ifdef TEST1
    {0.78E-01, 0.31E-01, 0.41E-02, 0.38E-02, 0.36E-02},
    {
        -0.450106022368E-02, -0.568426120793E-01, -0.876498238182E-01,
        0.531751834129E-01,  0.106350366826E+00,  0.531751834129E-01,
        0.192611490688E+00,  -0.331153469676E+00, -0.308498888708E+00,
        -0.616997777417E+00, -0.308498888708E+00, 0.645433404561E+00,
        -0.129198778272E+00, -0.258397556545E+00, -0.129198778272E+00,
        -0.119074508391E+01, -0.238149016782E+01, -0.119074508391E+01,
        -0.476298033563E+01, -0.238149016782E+01, -0.119074508391E+01,
    },
#else
#ifdef TEST2
    {0.17E+01, 0.17E+01, 0.81E-11, 0.81E-11, 0.81E-11},
    {
        -0.277344423214E+00, -0.902549628505E-01, -0.902549628505E-01,
        0.973880560793E-03,  0.194776112159E-02,  0.973880560793E-03,
        0.129726626281E-01,  -0.182704041173E-01, -0.381075641669E-03,
        -0.762151283338E-03, -0.381075641669E-03, 0.129726626281E-01,
        -0.381075641669E-03, -0.762151283338E-03, -0.381075641669E-03,
        -0.372455307820E-04, -0.744910615639E-04, -0.372455307820E-04,
        -0.148982123128E-03, -0.744910615639E-04, -0.372455307820E-04,
    }
#else
#ifdef TEST4
    {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06, 0.82E+06},
    {
        -0.258051271460E+01,
        -0.917654870253E-01,
        -0.924958459953E-01,
        0.504340587823E-06,
        0.100868117565E-05,
        0.504340587823E-06,
        /*
     v2rhoa2       = -0.673510067892E-03
     v2rhoab       = -0.140781632726E-02
     v2rhob2       = -0.681437355666E-03

     v2rhoasigmaaa =  0.880806595654E-08
     v2rhoasigmaab =  0.176161319131E-07
     v2rhoasigmabb =  0.880806595654E-08
     v2rhobsigmaaa =  0.892599758143E-08
     v2rhobsigmaab =  0.178519951629E-07
     v2rhobsigmabb =  0.892599758143E-08

     v2sigmaaa2    = -0.173276860156E-12
     v2sigmaaaab   = -0.346553720312E-12
     v2sigmaaabb   = -0.173276860156E-12
     v2sigmaab2    = -0.693107440625E-12
     v2sigmaabbb   = -0.346553720312E-12
     v2sigmabb2    = -0.173276860156E-12
        */
    }
#else
    {0.13E+00, 0.95E-01, 0.15E+00, 0.18E+00, 0.22E+00},
    {
        -0.475398805591E-02,
        -0.614930812180E-01,
        -0.682665101239E-01,
        0.486210103008E-02,
        0.972420206015E-02,
        0.486210103008E-02,
    }
/*
v2rhoa2       = -0.180150058142E+00
v2rhoab       = -0.388435750226E+00
v2rhob2       = -0.208827666674E+00

v2rhoasigmaaa =  0.288797891283E-01
v2rhoasigmaab =  0.577595782566E-01
v2rhoasigmabb =  0.288797891283E-01
v2rhobsigmaaa =  0.338842245007E-01
v2rhobsigmaab =  0.677684490015E-01
v2rhobsigmabb =  0.338842245007E-01

v2sigmaaa2    = -0.755929074438E-02
v2sigmaaaab   = -0.151185814888E-01
v2sigmaaabb   = -0.755929074438E-02
v2sigmaab2    = -0.302371629775E-01
v2sigmaabbb   = -0.151185814888E-01
v2sigmabb2    = -0.755929074438E-02
*/
#endif
#endif
#endif
};
