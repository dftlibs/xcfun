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

#pragma once

enum xc_functional_id {
  XC_SLATERX,
  XC_PW86X,
  XC_VWN3C,
  XC_VWN5C,
  XC_PBEC,
  XC_PBEX,
  XC_BECKEX,
  XC_BECKECORRX,
  XC_BECKESRX,
  XC_BECKECAMX,
  XC_BRX,
  XC_BRC,
  XC_BRXC,
  XC_LDAERFX,
  XC_LDAERFC,
  XC_LDAERFC_JT,
  XC_LYPC,
  XC_OPTX,
  XC_OPTXCORR,
  XC_REVPBEX,
  XC_RPBEX,
  XC_SPBEC,
  XC_VWN_PBEC,
  XC_KTX,
  XC_TFK,
  XC_TW,
  XC_PW91X,
  XC_PW91K,
  XC_PW92C,
  XC_M05X,
  XC_M05X2X,
  XC_M06X,
  XC_M06X2X,
  XC_M06LX,
  XC_M06HFX,
  XC_M05X2C,
  XC_M05C,
  XC_M06C,
  XC_M06HFC,
  XC_M06LC,
  XC_M06X2C,
  XC_TPSSC,
  XC_TPSSX,
  XC_REVTPSSC,
  XC_REVTPSSX,
  XC_PZ81C,
  XC_P86C,
  XC_P86CORRC,
  XC_BTK,
  XC_VWK,
  XC_B97X,
  XC_B97C,
  XC_B97_1X,
  XC_B97_1C,
  XC_B97_2X,
  XC_B97_2C,
  XC_CSC,
  XC_APBEC,
  XC_APBEX,
  XC_ZVPBESOLC,
  XC_BLOCX,
  XC_PBEINTC,
  XC_PBEINTX,
  XC_PBELOCC,
  XC_PBESOLX,
  XC_TPSSLOCC,
  XC_ZVPBEINTC,
  XC_PW91C,
  XC_NR_FUNCTIONALS
};

enum xc_parameter {
  XC_RANGESEP_MU = XC_NR_FUNCTIONALS,
  XC_EXX,
  XC_CAM_ALPHA,
  XC_CAM_BETA,
  XC_NR_PARAMETERS_AND_FUNCTIONALS
};
