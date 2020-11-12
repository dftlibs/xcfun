/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2020 Ulf Ekström and contributors.
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

namespace xcfun {
constexpr auto PROJECT_VERSION_MAJOR = 2;
constexpr auto PROJECT_VERSION_MINOR = 1;
constexpr auto PROJECT_VERSION_PATCH = 1;
constexpr auto XCFun_VERSION =
    ((PROJECT_VERSION_MAJOR << 16) | PROJECT_VERSION_MINOR | PROJECT_VERSION_PATCH);

auto version_as_string() noexcept -> std::string;

auto xcfun_get_version() noexcept -> unsigned int;
} // namespace xcfun
