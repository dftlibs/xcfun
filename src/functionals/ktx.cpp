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

template <typename num> static num ktx(const densvars<num> & d) {
  const parameter DELTA = 0.1;

  num ea = d.gaa / (DELTA + pow(d.a, 4.0 / 3.0));
  num eb = d.gbb / (DELTA + pow(d.b, 4.0 / 3.0));
  return ea + eb;
}

FUNCTIONAL(XC_KTX) = {
    "KT exchange GGA correction",
    "KT exchange GGA correction\n"
    "@article{keal:3015,\n"
    "author = {Thomas W. Keal and David J. Tozer},\n"
    "collaboration = {},\n"
    "title = {The exchange-correlation potential in Kohn--Sham nuclear magnetic "
    "resonance shielding calculations},\n"
    "publisher = {AIP},\n"
    "year = {2003},\n"
    "journal = {The Journal of Chemical Physics},\n"
    "volume = {119},\n"
    "number = {6},\n"
    "pages = {3015-3024},\n"
    "keywords = {eigenvalues and eigenfunctions; ab initio calculations; bond "
    "lengths; nuclear screening; nuclear magnetic resonance; density functional "
    "theory; ionisation potential; dissociation energies},\n"
    "url = {http://link.aip.org/link/?JCP/119/3015/1},\n"
    "doi = {10.1063/1.1590634}\n"
    "}\n"
    "xcfun version: Radovan Bast (radovan.bast@uit.no)\n"
    "tested against implementation in Dalton by Dave Wilson (davidwi@kjemi.uio.no)\n"
    "compared first derivatives only\n",
    XC_DENSITY | XC_GRADIENT,
    ENERGY_FUNCTION(ktx) XC_A_B_GAA_GAB_GBB,
    XC_PARTIAL_DERIVATIVES,
    2,
    1e-11,
    {0.39e+02, 0.38e+02, 0.81e+06, 0.82e+06, 0.82e+06},
    {
        0.12533312965365759737e+05, // 00000
        -.20906588852649329624e+03, // 10000
        -.22485950142470750279e+03, // 01000
        0.75553098007418890980e-02, // 00100
        0.00000000000000000000e+00, // 00001
        0.78213561302010112253e-02, // 00010
        0.12497415159317823097e+02, // 20000
        0.00000000000000000000e+00, // 11000
        -.25810603521789297204e-03, // 10100
        0.00000000000000000000e+00, // 10010
        0.00000000000000000000e+00, // 10001
        0.13794820570009042271e+02, // 02000
        0.00000000000000000000e+00, // 01100
        0.00000000000000000000e+00, // 01001
        -.27421890417647258121e-03, // 01010
        0.00000000000000000000e+00, // 00200
        0.00000000000000000000e+00, // 00110
        0.00000000000000000000e+00, // 00101
        0.00000000000000000000e+00, // 00020
        0.00000000000000000000e+00, // 00011
        0.00000000000000000000e+00  // 00002
    }};
