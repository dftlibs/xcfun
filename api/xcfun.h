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

#include <stdbool.h>
#include <stddef.h>

#include "XCFun/XCFunExport.h"

#define XCFun_API XCFun_EXPORT

#define XCFUN_API_VERSION 2

#ifdef __cplusplus
extern "C" {
#endif

#ifndef XC_MAX_ORDER
#define XC_MAX_ORDER 3
#endif

typedef enum {
  XC_MODE_UNSET = 0, // Need to be zero for default initialized structs
  XC_PARTIAL_DERIVATIVES,
  XC_POTENTIAL,
  XC_CONTRACTED,
  XC_NR_MODES
} xcfun_mode;

// Must be in sync with xcint_vars in xcint.cpp and with the fortran module
/*! \enum xcfun_vars
 *  \brief functional type
 */
// clang-format off
typedef enum {
  XC_VARS_UNSET = -1, /*!< Not defined */
  XC_A,               /*! LDA alpha */
  XC_N,               /*! LDA rho*/
  XC_A_B,             /*! LDA alpha & beta */
  XC_N_S,             /*! LDA rho and spin */

  XC_A_GAA,           /*! GGA with grad^2 alpha        */
  XC_N_GNN,           /*! GGA with grad^2 rho          */
  XC_A_B_GAA_GAB_GBB, /*! GGA with grad^2 alpha & beta */
  XC_N_S_GNN_GNS_GSS, /*! GGA with grad^2 rho and spin */
  XC_A_GAA_LAPA,      /*! metaGGA with grad^2 alpha        laplacian */
  XC_A_GAA_TAUA,      /*! metaGGA with grad^2 alpha        kinetic   */
  XC_N_GNN_LAPN,
  /*! metaGGA with grad^2 rho          laplacian */ // 10
  XC_N_GNN_TAUN,                /*! metaGGA with grad^2 rho          kinetic   */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB, /*! metaGGA with grad^2 alpha & beta laplacian */
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB, /*! metaGGA with grad^2 alpha & beta kinetic   */
  XC_N_S_GNN_GNS_GSS_LAPN_LAPS, /*! metaGGA with grad^2 rho and spin laplacian */
  XC_N_S_GNN_GNS_GSS_TAUN_TAUS, /*! metaGGA with grad^2 rho and spin kinetic   */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB, /*! metaGGA with grad^2 alpha & beta
                                             laplacian kinetic */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB, /*! metaGGA with grad^2 alpha &
                                                       beta laplacian kinetic current
                                                     */
  XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS, /*! metaGGA with grad^2 rho and spin
                                             laplacian kinetic */
  XC_A_AX_AY_AZ,            /*! GGA with gradient components alpha        */
  XC_A_B_AX_AY_AZ_BX_BY_BZ, /*! GGA with gradient components alpha & beta */
  XC_N_NX_NY_NZ,
  /*! GGA with gradient components rho          */ // 20
  XC_N_S_NX_NY_NZ_SX_SY_SZ, /*! GGA with gradient components rho and spin */
  XC_A_AX_AY_AZ_TAUA,       /*! metaGGA with gradient components alpha        */
  XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB, /*! metaGGA with gradient components alpha &
                                         beta */
  XC_N_NX_NY_NZ_TAUN, /*! metaGGA with gradient components rho          */
  XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS, /*! metaGGA with gradient components rho and
                                         spin */
  /* 2:nd order Taylor coefficients of alpha density, 1+3+6=10 numbers, rev gradlex
     order */
  XC_A_2ND_TAYLOR, /*! 2:nd order Taylor alpha        */
  /* 2:nd order Taylor expansion of alpha and beta densities (first alpha, then beta)
     20 numbers */
  XC_A_B_2ND_TAYLOR, /*! 2:nd order Taylor alpha & beta */
  XC_N_2ND_TAYLOR,   /*! 2:nd order Taylor rho          */
  XC_N_S_2ND_TAYLOR, /*! 2:nd order Taylor rho and spin */
  XC_NR_VARS
} xcfun_vars;
// clang-format on

XCFun_API const char * xcfun_version();

XCFun_API const char * xcfun_splash();

XCFun_API const char * xcfun_authors();

/*! \return The number of failed tests */
XCFun_API int xcfun_test();

/*! \brief Whether the library is compatible with the header file
 *  Checks that the compiled library and header file version match.
 *  Host should abort when that is not the case.
 *
 *  \warning This function should be called **before** instantiating
 *  any XCFunctional object.
 */
XCFun_API bool xcfun_is_compatible_library();

XCFun_API xcfun_vars xcfun_which_vars(const unsigned int func_type,
                                      const unsigned int dens_type,
                                      const unsigned int laplacian,
                                      const unsigned int kinetic,
                                      const unsigned int current,
                                      const unsigned int explicit_derivatives);

XCFun_API xcfun_mode xcfun_which_mode(const unsigned int mode_type);

XCFun_API const char * xcfun_enumerate_parameters(int param);

XCFun_API const char * xcfun_enumerate_aliases(int n);

XCFun_API const char * xcfun_describe_short(const char * name);

XCFun_API const char * xcfun_describe_long(const char * name);

/*! \struct xcfun_s
 *  Forward-declare opaque handle to a XCFunctional
 */
struct xcfun_s;

/*! \typedef xcfun_t
 *  Workaround to have xcfun_t available to C
 */
typedef struct xcfun_s xcfun_t;

XCFun_API xcfun_t * xcfun_new();

XCFun_API void xcfun_delete(xcfun_t * fun);

XCFun_API int xcfun_set(xcfun_t * fun, const char * name, double value);

XCFun_API int xcfun_get(const xcfun_t * fun, const char * name, double value[]);

XCFun_API bool xcfun_is_gga(const xcfun_t * fun);

XCFun_API bool xcfun_is_metagga(const xcfun_t * fun);

XCFun_API int xcfun_eval_setup(xcfun_t * fun,
                               xcfun_vars vars,
                               xcfun_mode mode,
                               int order);

XCFun_API int xcfun_user_eval_setup(xcfun_t * fun,
                                    const int order,
                                    const unsigned int func_type,
                                    const unsigned int dens_type,
                                    const unsigned int mode_type,
                                    const unsigned int laplacian,
                                    const unsigned int kinetic,
                                    const unsigned int current,
                                    const unsigned int explicit_derivatives);

XCFun_API int xcfun_input_length(const xcfun_t * fun);

XCFun_API int xcfun_output_length(const xcfun_t * fun);

XCFun_API void xcfun_eval(const xcfun_t * fun,
                          const double density[],
                          double result[]);

XCFun_API void xcfun_eval_vec(const xcfun_t * fun,
                              int nr_points,
                              const double * density,
                              int density_pitch,
                              double * result,
                              int result_pitch);
#ifdef __cplusplus
} // End of extern "C"
#endif
