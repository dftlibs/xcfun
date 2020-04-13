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

#include <array>

#include "XCFun/xcfun.h"
#include "functionals/list_of_functionals.hpp"

struct functional_data;

/*! \brief Exchange-correlation functional
 */
struct XCFunctional {
  explicit XCFunctional();

  int nr_active_functionals{0};
  int order{-1};
  int depends{0}; // XC_DENSITY, gradient etc
  xcfun_mode mode{XC_MODE_UNSET};
  xcfun_vars vars{XC_VARS_UNSET};
  std::array<functional_data *, XC_NR_FUNCTIONALS> active_functionals{{nullptr}};
  std::array<double, XC_NR_PARAMETERS_AND_FUNCTIONALS> settings{{0.0}};
};

namespace xcfun {
/*! Invalid order for given mode and vars */
constexpr auto XC_EORDER = 1;

/*! Invalid vars for functional type (ie. lda vars for gga) */
constexpr auto XC_EVARS = 2;

/*! Invalid mode for functional type (ie. potential for mgga) */
constexpr auto XC_EMODE = 4;

/// \cond DEV
XCFun_API XCFunctional * xcfun_new();
XCFun_API void xcfun_delete(XCFunctional *);
XCFun_API int xcfun_set(XCFunctional * fun, const char * name, double value);
XCFun_API int xcfun_get(const XCFunctional * fun, const char * name, double * value);
XCFun_API bool xcfun_is_gga(const XCFunctional * fun);
XCFun_API bool xcfun_is_metagga(const XCFunctional * fun);
XCFun_API int xcfun_eval_setup(XCFunctional * fun,
                               xcfun_vars vars,
                               xcfun_mode mode,
                               int order);
XCFun_API int xcfun_user_eval_setup(XCFunctional * fun,
                                    const int order,
                                    const unsigned int func_type,
                                    const unsigned int dens_type,
                                    const unsigned int mode_type,
                                    const unsigned int laplacian,
                                    const unsigned int kinetic,
                                    const unsigned int current,
                                    const unsigned int explicit_derivatives);
XCFun_API int xcfun_input_length(const XCFunctional * fun);
XCFun_API int xcfun_output_length(const XCFunctional * fun);
XCFun_API void xcfun_eval(const XCFunctional * fun,
                          const double input[],
                          double output[]);
XCFun_API void xcfun_eval_vec(const XCFunctional * fun,
                              int nr_points,
                              const double density[],
                              int density_pitch,
                              double result[],
                              int result_pitch);
/// \endcond
} // namespace xcfun
