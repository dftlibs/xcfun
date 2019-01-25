/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2018 Ulf Ekström and contributors.
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

#ifndef SWIG
#include "XCFun/XCFunExport.h"
#define XCFun_API XCFun_EXPORT
#else
#define XCFun_API
#endif

#define XCFUN_API_VERSION 2

#ifdef __cplusplus
extern "C" {
#endif

#ifndef XC_MAX_ORDER
#define XC_MAX_ORDER 3
#endif

// Used for regularizing input
#define XC_TINY_DENSITY 1e-14

#define XC_NO_REGULARIZATION

#define XC_EORDER 1 // Invalid order for given mode and vars
#define XC_EVARS 2  // Invalid vars for functional type (ie. lda vars for gga)
#define XC_EMODE 4  // Invalid mode for functional type (ie. potential for mgga)

enum xc_mode {
  XC_MODE_UNSET = 0, // Need to be zero for default initialized structs
  XC_PARTIAL_DERIVATIVES,
  XC_POTENTIAL,
  XC_CONTRACTED,
  XC_NR_MODES
};

// Must be in sync with xcint_vars in xcint.cpp and with the fortran module
/*! \enum xc_vars
 *  \brief functional type
 */
enum xc_vars {
  XC_VARS_UNSET = -1, /*!< Not defined */
  XC_A,               /*! LDA alpha */
  XC_N,               /*! LDA rho*/
  XC_A_B,             /*! LDA alpha & beta */
  XC_N_S,             /*! LDA rho and spin */

  XC_A_GAA,             /*! GGA with grad^2 alpha        */
  XC_N_GNN,             /*! GGA with grad^2 rho          */
  XC_A_B_GAA_GAB_GBB,   /*! GGA with grad^2 alpha & beta */
  XC_N_S_GNN_GNS_GSS,   /*! GGA with grad^2 rho and spin */
  XC_A_GAA_LAPA,                                    /*! metaGGA with grad^2 alpha        laplacian */
  XC_A_GAA_TAUA,                                    /*! metaGGA with grad^2 alpha        kinetic   */
  XC_N_GNN_LAPN,                                    /*! metaGGA with grad^2 rho          laplacian */ // 10
  XC_N_GNN_TAUN,                                    /*! metaGGA with grad^2 rho          kinetic   */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB,                     /*! metaGGA with grad^2 alpha & beta laplacian */
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,                     /*! metaGGA with grad^2 alpha & beta kinetic   */
  XC_N_S_GNN_GNS_GSS_LAPN_LAPS,                     /*! metaGGA with grad^2 rho and spin laplacian */
  XC_N_S_GNN_GNS_GSS_TAUN_TAUS,                     /*! metaGGA with grad^2 rho and spin kinetic   */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB,           /*! metaGGA with grad^2 alpha & beta laplacian kinetic */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB, /*! metaGGA with grad^2 alpha & beta laplacian kinetic current */
  XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS,           /*! metaGGA with grad^2 rho and spin laplacian kinetic */
  XC_A_AX_AY_AZ,             /*! GGA with gradient components alpha        */
  XC_A_B_AX_AY_AZ_BX_BY_BZ,  /*! GGA with gradient components alpha & beta */
  XC_N_NX_NY_NZ,             /*! GGA with gradient components rho          */ // 20
  XC_N_S_NX_NY_NZ_SX_SY_SZ,  /*! GGA with gradient components rho and spin */
  XC_A_AX_AY_AZ_TAUA,                 /*! metaGGA with gradient components alpha        */
  XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB, /*! metaGGA with gradient components alpha & beta */
  XC_N_NX_NY_NZ_TAUN,                 /*! metaGGA with gradient components rho          */
  XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS, /*! metaGGA with gradient components rho and spin */
  /* 2:nd order Taylor coefficients of alpha density, 1+3+6=10 numbers, rev gradlex order */
  XC_A_2ND_TAYLOR,    /*! 2:nd order Taylor alpha        */
  /* 2:nd order Taylor expansion of alpha and beta densities (first alpha, then beta) 20 numbers */
  XC_A_B_2ND_TAYLOR,  /*! 2:nd order Taylor alpha & beta */
  XC_N_2ND_TAYLOR,    /*! 2:nd order Taylor rho          */
  XC_N_S_2ND_TAYLOR,  /*! 2:nd order Taylor rho and spin */
  XC_NR_VARS
};

XCFun_API double xcfun_version(void);
XCFun_API const char * xcfun_splash(void);
XCFun_API int xcfun_test(void);

typedef struct xc_functional_obj * xc_functional;

XCFun_API xc_functional xc_new_functional_not_macro(int api_version);
#define xc_new_functional() xc_new_functional_not_macro(XCFUN_API_VERSION)

XCFun_API void xc_free_functional(xc_functional fun);

// Fill the data array with information to recreate an exact copy of fun.
// This is not stable between different xcfun "versions"
// If buflen > 0 return the number of elements written,
// else return the number of elements needed. buf is not accessed if buflen <= 0.
int xc_serialize(xc_functional fun, int buflen, double * buf);
// make fun an exact copy of the functional used to fill buf (with xc_serialize)
void xc_deserialize(xc_functional fun, double * buf);

// Call with n >= 0, Returns a pointer to an internal string with the name
// of parameter nr n. Return NULL when n is too large.
XCFun_API const char * xc_enumerate_parameters(int n);
// Like xc_enumerate_parameters, but over aliases
XCFun_API const char * xc_enumerate_aliases(int n);
// Try to either set or get a parameter. Return 0 if all was well,
// otherwise the name was invalid.
XCFun_API int xc_set(xc_functional fun, const char * name, double value);
XCFun_API int xc_get(xc_functional fun, const char * name, double * value);
XCFun_API const char * xc_describe_short(const char * name);
const char * xc_describe_long(const char * name);

XCFun_API int xc_is_gga(xc_functional fun);
XCFun_API int xc_is_metagga(xc_functional fun);

XCFun_API int xc_set_fromstring(xc_functional fun, const char * str); // Defines a functional
                                                            // from a string on the
                                                            // form "fun[=value]"

XCFun_API int xc_user_eval_setup(xc_functional fun,
                       const int order, // order of the derivative requested (order=1 is the xc potential)
                       const unsigned int func_type, // LDA (0), GGA (1), metaGGA (2), taylor (3)
                       const unsigned int dens_type,  // A (0), N (1), A_B (2), N_S (3)
                       const unsigned int mode_type,  // same as the enum list
                       const unsigned int laplacian,  // 0/1 laplacian no/yes
                       const unsigned int kinetic,    // 0/1 kinetic energy no/yes
                       const unsigned int current,    // 0/1 current density no/yes
                       const unsigned int explicit_derivatives);   // 0/1 gamma vs explicit partial derivatives

// Try to set the functional evaluation vars, mode and order
// return some combination of XC_E* if an error occurs, else 0.
XCFun_API int xc_eval_setup(xc_functional fun,
                            enum xc_vars vars,
                            enum xc_mode mode,
                            int order);

// Length of the density[] argument to eval()
XCFun_API int xc_input_length(xc_functional fun);

// Length of the result[] argument to eval()
XCFun_API int xc_output_length(xc_functional fun);

// Evaluate the functional at density
// In contracted mode density is of dimension 2^porder*Nvars
XCFun_API void xc_eval(xc_functional fun, const double * density, double * result);
/* Vector version of xc_eval.
   density_pitch = density[start_of_second_point] - density[start_of_first_point],
   likewise for result_pitch. */
XCFun_API void xc_eval_vec(xc_functional fun,
                 int nr_points,
                 const double * density,
                 int density_pitch,
                 double * result,
                 int result_pitch);

// Index into result[] for derivative with given index (length as input_length() )
int xc_derivative_index(xc_functional fun, const int derivative[]);

#ifdef __cplusplus
} // End of extern "C"
#endif

// Derivative indices into xc_eval output in partial derivative mode

#define XC_D0 0
#define XC_D1 1
#define XC_D2 2

#define XC_D00 0
#define XC_D10 1
#define XC_D01 2
#define XC_D20 3
#define XC_D11 4
#define XC_D02 5

#define XC_D00000 0
#define XC_D10000 1
#define XC_D01000 2
#define XC_D00100 3
#define XC_D00010 4
#define XC_D00001 5
#define XC_D20000 6
#define XC_D11000 7
#define XC_D10100 8
#define XC_D10010 9
#define XC_D10001 10
#define XC_D02000 11
#define XC_D01100 12
#define XC_D01010 13
#define XC_D01001 14
#define XC_D00200 15
#define XC_D00110 16
#define XC_D00101 17
#define XC_D00020 18
#define XC_D00011 19
#define XC_D00002 20
