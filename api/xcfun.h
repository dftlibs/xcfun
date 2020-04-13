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

#define XCFun_API XCFUN_EXPORT

#define XCFUN_API_VERSION 2 /*!< Version of the XCFun API */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef XCFUN_MAX_ORDER
#define XCFUN_MAX_ORDER 6 /*!< Maximum differentiation order for XC kernels */
#endif

/*! \brief Evaluation mode for functional derivatives */
typedef enum {
  XC_MODE_UNSET = 0,      /*!< Need to be zero for default initialized structs */
  XC_PARTIAL_DERIVATIVES, /*!< ??? */
  XC_POTENTIAL,           /*!< ??? */
  XC_CONTRACTED,          /*!< ??? */
  XC_NR_MODES             /*!< ??? */
} xcfun_mode;

// Must be in sync with xcint_vars in xcint.cpp and with the fortran module
// clang-format off
/*! \brief Types of variables to define a functional.
 *
 *  The XC energy density and derivatives can be evaluated using a variety of
 *  variables and variables combinations. The variables in this `enum` are named as:
 *
 *    - `XC_` prefix
 *    - Tag for density variables.
 *    - Tag for gradient variables.
 *    - Tag for Laplacian variables.
 *    - Tag for kinetic energy density variables.
 *    - Tag for current density variables.
 *
 * XCFun recognizes the following basic variables:
 *
 *    - `A`, the spin-up electron number density: \f$n_{\alpha}\f$
 *    - `B`, the spin-down electron number density: \f$n_{\beta}\f$
 *    - `GAA`, the square magnitude of the spin-up density gradient: \f$\sigma_{\alpha \alpha} = \nabla n_\alpha.\nabla n_\alpha\f$
 *    - `GAB`, the dot product of the spin-up and spin-down density gradients: \f$\sigma_{\alpha \beta} = \nabla n_\alpha.\nabla n_\beta\f$
 *    - `GBB`, the square magnitude of the spin-down density gradient: \f$\sigma_{\beta \beta} = \nabla n_\beta.\nabla n_\beta\f$
 *    - `LAPA`, the Laplacian of the spin-up density: \f$\nabla^2 n_{\alpha}\f$
 *    - `LAPB`, the Laplacian of the spin-down density: \f$\nabla^2 n_{\beta}\f$
 *    - `TAUA`, the spin-up Kohn-Sham kinetic energy density: \f$\tau_\alpha = \frac{1}{2} \sum_i |\psi_{i \alpha}|^2\f$
 *    - `TAUB`, the spin-down Kohn-Sham kinetic energy density: \f$\tau_\beta = \frac{1}{2} \sum_i |\psi_{i \beta}|^2\f$
 *    - `JPAA`, the spin-up current density: \f$\mathbf{j}_{\alpha\alpha}\f$
 *    - `JPBB`, the spin-down current density: \f$\mathbf{j}_{\beta\beta}\f$
 *
 *  The following quantities are also recognized:
 *
 *    - `N`, the number density: \f$n = n_{\alpha} + n_{\beta}\f$
 *    - `S`, the spin density: \f$s = n_{\alpha} - n_{\beta}\f$
 *    - `GNN`, the square magnitude of the density gradient: \f$\sigma_{nn} = \nabla n.\nabla n\f$
 *    - `GSS`, the dot product of the number and spin density gradients: \f$\sigma_{ns} = \nabla n.\nabla s\f$
 *    - `GNS`, the square magnitude of the spin density gradient: \f$\sigma_{ss} = \nabla s.\nabla s\f$
 *    - `LAPN`, the Laplacian of the density: \f$\nabla^2 n\f$
 *    - `LAPS`, the Laplacian of the spin density: \f$\nabla^2 s\f$
 *    - `TAUN`, the Kohn-Sham kinetic energy density: \f$\tau_n\f$
 *    - `TAUS`, the spin Kohn-Sham kinetic energy density: \f$\tau_s\f$
 *
 *  XC functionals depending on the gradient of the density can furthermore be
 *  defined to use the \f$(x, y, z)\f$ components of the gradient explicitly.
 */
typedef enum {
  XC_VARS_UNSET = -1, /*!< Not defined */
  XC_A,               /*!< LDA with \f$n_{\alpha}\f$ */
  XC_N,               /*!< LDA with \f$n\f$ */
  XC_A_B,             /*!< LDA with \f$n_{\alpha}\f$ and \f$n_{\beta}\f$ */
  XC_N_S,             /*!< LDA with \f$n\f$ and \f$s\f$ */

  XC_A_GAA,             /*!< GGA with grad^2 alpha        */
  XC_N_GNN,             /*!< GGA with grad^2 rho          */
  XC_A_B_GAA_GAB_GBB,   /*!< GGA with grad^2 alpha & beta */
  XC_N_S_GNN_GNS_GSS,   /*!< GGA with grad^2 rho and spin */
  XC_A_GAA_LAPA,                                    /*!< metaGGA with grad^2 alpha        laplacian */
  XC_A_GAA_TAUA,                                    /*!< metaGGA with grad^2 alpha        kinetic   */
  XC_N_GNN_LAPN,                                    /*!< metaGGA with grad^2 rho          laplacian */ // 10
  XC_N_GNN_TAUN,                                    /*!< metaGGA with grad^2 rho          kinetic   */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB,                     /*!< metaGGA with grad^2 alpha & beta laplacian */
  XC_A_B_GAA_GAB_GBB_TAUA_TAUB,                     /*!< metaGGA with grad^2 alpha & beta kinetic   */
  XC_N_S_GNN_GNS_GSS_LAPN_LAPS,                     /*!< metaGGA with grad^2 rho and spin laplacian */
  XC_N_S_GNN_GNS_GSS_TAUN_TAUS,                     /*!< metaGGA with grad^2 rho and spin kinetic   */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB,           /*!< metaGGA with grad^2 alpha & beta laplacian kinetic */
  XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB, /*!< metaGGA with grad^2 alpha & beta laplacian kinetic current */
  XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS,           /*!< metaGGA with grad^2 rho and spin laplacian kinetic */
  XC_A_AX_AY_AZ,             /*!< GGA with gradient components alpha        */
  XC_A_B_AX_AY_AZ_BX_BY_BZ,  /*!< GGA with gradient components alpha & beta */
  XC_N_NX_NY_NZ,             /*!< GGA with gradient components rho          */ // 20
  XC_N_S_NX_NY_NZ_SX_SY_SZ,  /*!< GGA with gradient components rho and spin */
  XC_A_AX_AY_AZ_TAUA,                 /*!< metaGGA with gradient components alpha        */
  XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB, /*!< metaGGA with gradient components alpha & beta */
  XC_N_NX_NY_NZ_TAUN,                 /*!< metaGGA with gradient components rho          */
  XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS, /*!< metaGGA with gradient components rho and spin */

  XC_A_2ND_TAYLOR,    /*!< 2nd order Taylor coefficients of alpha density, 1+3+6=10 numbers, rev gradlex order */
  XC_A_B_2ND_TAYLOR,  /*!< 2nd order Taylor expansion of alpha and beta densities (first alpha, then beta) 20 numbers */
  XC_N_2ND_TAYLOR,    /*!< 2nd order Taylor rho          */
  XC_N_S_2ND_TAYLOR,  /*!< 2nd order Taylor rho and spin */
  XC_NR_VARS          /*!< Number of variables */
} xcfun_vars;
// clang-format on

/*! \brief The version of XCFun in use
 *  \return the version of XCFun
 */
XCFun_API const char * xcfun_version();

/*! \brief The XCFun splash screen
 *  \return A `char` array with the XCFun splash screen.
 *
 *  Return a multi-line string describing the library. This functions shows the
 *  code attribution and literature citation.
 *  It should be called when initializing XCFun in client code, so that your
 *  users find the right citation for the library.
 */
XCFun_API const char * xcfun_splash();

/*! \brief The XCFun splash screen
 *  \return A `char` array with the current list of XCFun authors.
 */
XCFun_API const char * xcfun_authors();

/*! \brief Test XCFun
 *  \return the number of failed tests.
 *
 *  Run all internal tests and return the number of failed tests.
 */
XCFun_API int xcfun_test();

/*! \brief Whether the library is compatible with the header file
 *  Checks that the compiled library and header file version match.
 *  Host should abort when that is not the case.
 *
 *  \warning This function should be called **before** instantiating
 *  any XCFunctional object.
 */
XCFun_API bool xcfun_is_compatible_library();

// clang-format off
/*! \brief Obtain correct value of `xcfun_vars` `enum`.
 *  \param[in] func_type LDA (0), GGA (1), metaGGA (2), taylor (3)
 *  \param[in] dens_type Alpha (A,0), Rho (N,1), Alpha&Beta (A_B,2), Rho&Spin (N_S,3)
 *  \param[in] laplacian (0 not required / 1 required)
 *  \param[in] kinetic  (0 not required / 1 required)
 *  \param[in] current   (0 not required / 1 required)
 *  \param[in] explicit_derivatives  (0 not required / 1 required)
 *  \return XC functional variables to use
 *
 *  This routine encodes the different options bitwise. Each legitimate
 *  combination is then converted to the corresponding enum value.
 *
 *  \rst
 *
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 |                                                |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  | 0 | 0 |   |   |   |   |   |   | LDA                                            |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  | 0 | 1 |   |   |   |   |   |   | GGA                                            |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  | 1 | 0 |   |   |   |   |   |   | metaGGA                                        |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  | 1 | 1 |   |   |   |   |   |   | Taylor                                         |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   | 0 | 0 |   |   |   |   | :math:`\rho_{\alpha}`                          |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   | 0 | 1 |   |   |   |   | :math:`\rho`                                   |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   | 1 | 0 |   |   |   |   | :math:`\rho_{\alpha}` and :math:`\rho_{\beta}` |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   | 1 | 1 |   |   |   |   | :math:`\rho` and :math:`s`                     |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   | 0 |   |   |   | no laplacian                                   |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   | 1 |   |   |   | laplacian required                             |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   |   | 0 |   |   | no kinetic energy density                      |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   |   | 1 |   |   | kinetic energy density required                |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   |   |   | 0 |   | no current density required                    |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   |   |   | 1 |   | current density required                       |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   |   |   |   | 0 | :math:`\gamma`-type partial derivatives        |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *  |   |   |   |   |   |   |   | 1 | explicit partial derivatives                   |
 *  +---+---+---+---+---+---+---+---+------------------------------------------------+
 *
 *  \endrst
 */
// clang-format on
XCFun_API xcfun_vars xcfun_which_vars(const unsigned int func_type,
                                      const unsigned int dens_type,
                                      const unsigned int laplacian,
                                      const unsigned int kinetic,
                                      const unsigned int current,
                                      const unsigned int explicit_derivatives);

/*! \brief Obtain correct value of `xcfun_mode` `enum`.
 *  \param[in] mode_type Partial derivatives (1), Potential (2), Contracted (3)
 *  \return The XC functional evaluation mode
 */
XCFun_API xcfun_mode xcfun_which_mode(const unsigned int mode_type);

/*! \brief Describe XC functional parameters
 *  \param[in] param the parameter to describe. `param` >= 0.
 *  \return description of the given parameter, or `NULL` is `param` is too large.
 */
XCFun_API const char * xcfun_enumerate_parameters(int param);

/*! \brief Describe XC functional aliases
 *  \param[in] n the alias to describe. `n` >= 0.
 *  \return description of the given alias, or `NULL` is `n` is too large.
 */
XCFun_API const char * xcfun_enumerate_aliases(int n);

/*! \brief Short description of the XC functional
 *  \param[in] name
 *  \return short description of the functional.
 */
XCFun_API const char * xcfun_describe_short(const char * name);

/*! \brief Long description of the XC functional
 *  \param[in] name
 *  \return long description of the functional.
 */
XCFun_API const char * xcfun_describe_long(const char * name);

/*! \struct xcfun_s
 *  Forward-declare opaque handle to a `XCFunctional` object.
 */
struct xcfun_s;

/*! \typedef xcfun_t
 *  \brief Opaque handle to a `XCFunctional` object.
 *  \note This type definition is a workaround to have the opaque `xcfun_t`
 *  `struct` available to C.
 */
typedef struct xcfun_s xcfun_t;

/*! \brief Create a new XC functional object
 *  \return A `xcfun_t` object.
 *
 *  Create a new functional object. The creation of this
 *  object may be rather slow; create an object once for each calculation, not
 *  once for each grid point.
 */
XCFun_API xcfun_t * xcfun_new();

/*! \brief Delete a XCFun functional
 *  \param[in, out] fun the XCFun functional to be deleted
 */
XCFun_API void xcfun_delete(xcfun_t * fun);

/*! \brief Set a parameter in the XC functional
 *  \param[in, out] fun
 *  \param[in] name
 *  \param[in] value
 *  \return error code (0 means normal exit)
 */
XCFun_API int xcfun_set(xcfun_t * fun, const char * name, double value);

/*! \brief Get weight of given functional in the current setup
 *  \param[in] fun the functional object
 *  \param[in] name functional name to test, aliases not supported
 *  \param[out] value weight of functional
 *
 *  \return `0` if `name` is a valid functional, `-1` if not.
 *  See `list_of_functionals.hpp` for valid functional names.
 */
XCFun_API int xcfun_get(const xcfun_t * fun, const char * name, double * value);

/*! \brief Is the XC functional GGA?
 *  \param[in, out] fun
 *  \return Whether `fun` is a GGA-type functional
 */
XCFun_API bool xcfun_is_gga(const xcfun_t * fun);

/*! \brief Is the XC functional GGA?
 *  \param[in, out] fun
 *  \return Whether `fun` is a metaGGA-type functional
 */
XCFun_API bool xcfun_is_metagga(const xcfun_t * fun);

/*! \brief Set up XC functional evaluation variables, mode, and order
 *  \param[in, out] fun XC functional object
 *  \param[in] vars evaluation variables
 *  \param[in] mode evaluation mode
 *  \param[in] order order of the derivative requested (order=1 is the xc potential)
 *  \return some combination of `XC_E*` if an error occurs, else 0
 */
XCFun_API int xcfun_eval_setup(xcfun_t * fun,
                               xcfun_vars vars,
                               xcfun_mode mode,
                               int order);

/*! \brief Host program-friendly set up of the XC functional evaluation variables,
 * mode, and order
 *  \param[in, out] fun XC functional object \param[in] order order of the derivative
 requested (order 0 (functional), 1 (potential), 2 (hessian), ....)
 *  \param[in] func_type LDA (0), GGA (1), metaGGA (2), taylor (3)
 *  \param[in] dens_type Alpha (A,0), Rho (N,1), Alpha&Beta (A_B,2), Rho&Spin (N_S,3)
 *  \param[in] mode_type Partial derivatives (1), Potential (2), Contracted (3)
 *  \param[in] laplacian (0 not required / 1 required)
 *  \param[in] kinetic  (0 not required / 1 required)
 *  \param[in] current   (0 not required / 1 required)
 *  \param[in] explicit_derivatives  (0 not required / 1 required)
 *  \return some combination of `XC_E*` if an error occurs, else 0
 */
XCFun_API int xcfun_user_eval_setup(xcfun_t * fun,
                                    const int order,
                                    const unsigned int func_type,
                                    const unsigned int dens_type,
                                    const unsigned int mode_type,
                                    const unsigned int laplacian,
                                    const unsigned int kinetic,
                                    const unsigned int current,
                                    const unsigned int explicit_derivatives);

/*! \brief Length of the density[] argument to `xcfun_eval`
 *  \param[in, out] fun XC functional object
 *  \return some combination of `XC_E*` if an error occurs, else 0
 */
XCFun_API int xcfun_input_length(const xcfun_t * fun);

/*! \brief Length of the result[] argument to `xcfun_eval`
 *  \param[in, out] fun XC functional object
 *  \return Return the number of output coefficients computed by `xc_eval()`.
 *
 *  \note All derivatives up to order are calculated, not only those of the
 * particular order.
 */
XCFun_API int xcfun_output_length(const xcfun_t * fun);

/*! \brief Evaluate the XC functional for given density at a point.
 *  \param[in, out] fun XC functional object
 *  \param[in] density
 *  \param[in, out] result
 *
 *  \note In contracted mode density is of dimension
 * \f$2^{\mathrm{order}}*N_{\mathrm{vars}}\f$
 */
XCFun_API void xcfun_eval(const xcfun_t * fun,
                          const double density[],
                          double result[]);

/*! \brief Evaluate the XC functional for given density on a set of points.
 *  \param[in, out] fun XC functional object
 *  \param[in] nr_points number of points in the evaluation set.
 *  \param[in] density
 *  \param[in] density_pitch `density[start_of_second_point] -
 * density[start_of_first_point]` \param[in, out] result
 *  \param[in] result_pitch
 * `result[start_of_second_point] - result[start_of_first_point]`
 *
 *  \note In contracted mode density is of dimension
 * \f$2^{\mathrm{order}}*N_{\mathrm{vars}}\f$
 */
XCFun_API void xcfun_eval_vec(const xcfun_t * fun,
                              int nr_points,
                              const double * density,
                              int density_pitch,
                              double * result,
                              int result_pitch);
#ifdef __cplusplus
} // End of extern "C"
#endif
