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

PARAMETER(XC_RANGESEP_MU) = {"Range separation inverse length [1/a0]", 0.4};

PARAMETER(XC_EXX) = {
    "Amount of exact (HF like) exchange (must be provided externally)",
    0.0};

PARAMETER(XC_CAM_ALPHA) = {
    "Amount of exact (HF like) exchange within CAM-B3LYP functional",
    0.19};

PARAMETER(XC_CAM_BETA) = {
    "Amount of long-range (HF like) exchange within CAM-B3LYP functional",
    0.46};
