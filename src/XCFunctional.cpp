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

#include "XCFunctional.hpp"
#include "XCFun/xcfun.h"

#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>

#include "functionals/list_of_functionals.hpp"
#include "version_info.hpp"
#include "xcint.hpp"

#ifndef AS_TYPE
#define AS_TYPE(Type, Obj) reinterpret_cast<Type *>(Obj)
#endif

#ifndef AS_CTYPE
#define AS_CTYPE(Type, Obj) reinterpret_cast<const Type *>(Obj)
#endif

namespace xcfun {
/*! Invalid order for given mode and vars */
constexpr auto XC_EORDER = 1;

/*! Invalid vars for functional type (ie. lda vars for gga) */
constexpr auto XC_EVARS = 2;

/*! Invalid mode for functional type (ie. potential for mgga) */
constexpr auto XC_EMODE = 4;

auto version_as_string() noexcept -> std::string {
  std::ostringstream stream;
  stream << std::to_string(PROJECT_VERSION_MAJOR) << "."
         << std::to_string(PROJECT_VERSION_MINOR) << "."
         << std::to_string(PROJECT_VERSION_PATCH);
  return stream.str();
}

auto xcfun_get_version() noexcept -> unsigned int { return XCFun_VERSION; }
} // namespace xcfun

const char * xcfun_version() { return xcfun::version_as_string().c_str(); }

const char * xcfun_splash() {
  return "XCFun DFT library Copyright 2009-2020 Ulf Ekstrom and contributors.\n"
         "See http://dftlibs.org/xcfun/ for more information.\n\n"
         "This is free software; see the source code for copying conditions.\n"
         "There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or\n"
         "FITNESS FOR A PARTICULAR PURPOSE. For details see the documentation.\n"
         "Scientific users of this library should cite\n"
         "U. Ekstrom, L. Visscher, R. Bast, A. J. Thorvaldsen and K. Ruud;\n"
         "J.Chem.Theor.Comp. 2010, DOI: 10.1021/ct100117s\n";
}

const char * xcfun_authors() {
  return "XCFun was written by Ulf Ekstrom, with contributions from\n"
         "Andre S. P. Gomes\n"
         "Radovan Bast\n"
         "Andrea Debnarova\n"
         "Paola Gori-Giorgi\n"
         "Alexei Yakovlev\n"
         "Michael Seth\n";
}

int xcfun_test() {
  int nfail = 0, res;
  for (auto f = 0; f < XC_NR_FUNCTIONALS; ++f) {
    const functional_data * fd = &xcint_funs[f];
    auto fun = xcfun_new();
    xcfun_set(fun, fd->name, 1.0);

    if (fd->test_mode != XC_MODE_UNSET) {
      if ((res = xcfun_eval_setup(
               fun, fd->test_vars, fd->test_mode, fd->test_order)) == 0) {
        int n = xcfun_output_length(fun);
        auto out = new double[n];
        if (!fd->test_in)
          xcfun::die("Functional has no test input!", f);
        xcfun_eval(fun, fd->test_in, out);
        int nerr = 0;
        for (auto i = 0; i < n; ++i)
          if (std::abs(out[i] - fd->test_out[i]) >
              std::abs(fd->test_out[i] * fd->test_threshold))
            nerr++;
        if (nerr > 0) {
          fprintf(stderr,
                  "Error detected in functional %s with tolerance %g:\n",
                  fd->name,
                  fd->test_threshold);
          fprintf(stderr, "Abs.Error \tComputed              Reference\n");
          for (auto i = 0; i < n; ++i) {
            fprintf(stderr, "%.1e", std::abs(out[i] - fd->test_out[i]));
            fprintf(stderr, "    %+.16e \t%+.16e", out[i], fd->test_out[i]);
            if (std::abs(out[i] - fd->test_out[i]) >
                std::abs(fd->test_out[i] * fd->test_threshold))
              fprintf(stderr, " *");
            fprintf(stderr, "\n");
          }
          nfail++;
        }
        delete[] out;
      } else {
        fprintf(stderr,
                "Functional %s not supporting its own test, error %i\n",
                fd->name,
                res);
        nfail++;
      }
    } else {
      fprintf(stderr, "%s has no test\n", fd->name);
    }
    xcfun_delete(fun);
  }
  return nfail;
}

bool xcfun_is_compatible_library() {
  unsigned int major = (xcfun::xcfun_get_version() >> 16);
  return (major == xcfun::PROJECT_VERSION_MAJOR);
}

xcfun_vars xcfun_which_vars(const unsigned int func_type,
                            const unsigned int dens_type,
                            const unsigned int laplacian,
                            const unsigned int kinetic,
                            const unsigned int current,
                            const unsigned int explicit_derivatives) {
  if (func_type > 3 || dens_type > 3 || laplacian > 1 || kinetic > 1 ||
      current > 1 || explicit_derivatives > 1) {
    xcfun::die("xcfun_which_vars: invalid input", -1);
  }

  xcfun_vars vars = XC_VARS_UNSET;

  /* Bit encoding of information
     7   6   5   4   3   2   1   0
  ------------------------------------------------------------------
     0   0                             LDA
     0   1                             GGA
     1   0                             metaGGA
     1   1                             Taylor
  ------------------------------------------------------------------
             0   0                     alpha density
             0   1                     n density
             1   0                     anpha and beta densities
             1   1                     n and s density
  ------------------------------------------------------------------
                     0                 no laplacian
                     1                 laplacian required
  ------------------------------------------------------------------
                         0             no kinetic energy
                         1             kinetic energy required
  ------------------------------------------------------------------
                             0         no current density required
                             1         current density required
  ------------------------------------------------------------------
                                 0     gamma-type partial derivatives
                                 1     explicit partial derivatives
   */

  int bitwise_vars = 0;
  bitwise_vars += func_type << 6;       // 0 64 128 192
  bitwise_vars += dens_type << 4;       // 0 16  32  48
  bitwise_vars += laplacian << 3;       // 0  8
  bitwise_vars += kinetic << 2;         // 0  4
  bitwise_vars += current << 1;         // 0  2
  bitwise_vars += explicit_derivatives; // 0  1

  switch (bitwise_vars) {
    case (0):
      vars = XC_A;
      break; // 0  0  |  0  0  |  0  |  0  |  0  |  0
    case (16):
      vars = XC_N;
      break; // 0  0  |  0  1  |  0  |  0  |  0  |  0
    case (32):
      vars = XC_A_B;
      break; // 0  0  |  1  0  |  0  |  0  |  0  |  0
    case (48):
      vars = XC_N_S;
      break; // 0  0  |  1  1  |  0  |  0  |  0  |  0
    case (64):
      vars = XC_A_GAA;
      break; // 0  1  |  0  0  |  0  |  0  |  0  |  0
    case (65):
      vars = XC_A_AX_AY_AZ;
      break; // 0  1  |  0  0  |  0  |  0  |  0  |  1
    case (80):
      vars = XC_N_GNN;
      break; // 0  1  |  0  1  |  0  |  0  |  0  |  0
    case (81):
      vars = XC_N_NX_NY_NZ;
      break; // 0  1  |  0  1  |  0  |  0  |  0  |  1
    case (96):
      vars = XC_A_B_GAA_GAB_GBB;
      break; // 0  1  |  1  0  |  0  |  0  |  0  |  0
    case (97):
      vars = XC_A_B_AX_AY_AZ_BX_BY_BZ;
      break; // 0  1  |  1  0  |  0  |  0  |  0  |  1
    case (112):
      vars = XC_N_S_GNN_GNS_GSS;
      break; // 0  1  |  1  1  |  0  |  0  |  0  |  0
    case (113):
      vars = XC_N_S_NX_NY_NZ_SX_SY_SZ;
      break; // 0  1  |  1  1  |  0  |  0  |  0  |  1
    case (132):
      vars = XC_A_GAA_TAUA;
      break; // 1  0  |  0  0  |  0  |  1  |  0  |  0
    case (133):
      vars = XC_A_AX_AY_AZ_TAUA;
      break; // 1  0  |  0  0  |  0  |  1  |  0  |  1
    case (136):
      vars = XC_A_GAA_LAPA;
      break; // 1  0  |  0  0  |  1  |  0  |  0  |  0
    case (148):
      vars = XC_N_GNN_TAUN;
      break; // 1  0  |  0  1  |  0  |  1  |  0  |  0
    case (149):
      vars = XC_N_NX_NY_NZ_TAUN;
      break; // 1  0  |  0  1  |  0  |  1  |  0  |  1
    case (152):
      vars = XC_N_GNN_LAPN;
      break; // 1  0  |  0  1  |  1  |  0  |  0  |  0
    case (164):
      vars = XC_A_B_GAA_GAB_GBB_TAUA_TAUB;
      break; // 1  0  |  1  0  |  0  |  1  |  0  |  0
    case (165):
      vars = XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB;
      break; // 1  0  |  1  0  |  0  |  1  |  0  |  1
    case (168):
      vars = XC_A_B_GAA_GAB_GBB_LAPA_LAPB;
      break; // 1  0  |  1  0  |  1  |  0  |  0  |  0
    case (172):
      vars = XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB;
      break; // 1  0  |  1  0  |  1  |  1  |  0  |  0
    case (174):
      vars = XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB;
      break; // 1  0  |  1  0  |  1  |  1  |  1  |  0
    case (180):
      vars = XC_N_S_GNN_GNS_GSS_TAUN_TAUS;
      break; // 1  0  |  1  1  |  0  |  1  |  0  |  0
    case (181):
      vars = XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS;
      break; // 1  0  |  1  1  |  0  |  1  |  0  |  1
    case (184):
      vars = XC_N_S_GNN_GNS_GSS_LAPN_LAPS;
      break; // 1  0  |  1  1  |  1  |  0  |  0  |  0
    case (188):
      vars = XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS;
      break; // 1  0  |  1  1  |  1  |  1  |  0  |  0
    case (192):
      vars = XC_A_2ND_TAYLOR;
      break; // 1  1  |  0  0  |  0  |  0  |  0  |  0
    case (208):
      vars = XC_N_2ND_TAYLOR;
      break; // 1  1  |  0  1  |  0  |  0  |  0  |  0
    case (224):
      vars = XC_A_B_2ND_TAYLOR;
      break; // 1  1  |  1  0  |  0  |  0  |  0  |  0
    case (240):
      vars = XC_N_S_2ND_TAYLOR;
      break; // 1  1  |  1  1  |  0  |  0  |  0  |  0
    default:
      xcfun::die("xc_user_eval_setup: Invalid vars", bitwise_vars);
  }

  return vars;
}

xcfun_mode xcfun_which_mode(const unsigned int mode_type) {
  if (mode_type > 3) {
    xcfun::die("xcfun_which_mode: invalid input", -1);
  }
  xcfun_mode mode = XC_MODE_UNSET;

  switch (mode_type) {
    case (1):
      mode = XC_PARTIAL_DERIVATIVES;
      break;
    case (2):
      mode = XC_POTENTIAL;
      break;
    case (3):
      mode = XC_CONTRACTED;
      break;
    default:
      xcfun::die("xc_user_eval_setup: Invalid mode", mode_type);
  }

  return mode;
}

const char * xcfun_enumerate_parameters(int param) {
  xcint_assure_setup();
  if (param >= 0 && param < XC_NR_FUNCTIONALS)
    return xcint_funs[param].name;
  else if (param < XC_NR_PARAMETERS_AND_FUNCTIONALS)
    return xcint_params[param].name;
  else
    return 0;
}

// Like xcfun_enumerate_parameters, but over aliases
const char * xcfun_enumerate_aliases(int n) {
  if (n >= 0 and n < XC_MAX_ALIASES)
    return xcint_aliases[n].name;
  else
    return 0;
}

const char * xcfun_describe_short(const char * name) {
  xcint_assure_setup();
  int k;
  if ((k = xcint_lookup_functional(name)) >= 0)
    return xcint_funs[k].short_description;
  else if ((k = xcint_lookup_parameter(name)) >= 0)
    return xcint_params[k].description;
  else if ((k = xcint_lookup_alias(name)) >= 0)
    return xcint_aliases[k].description;
  else
    return 0;
}

const char * xcfun_describe_long(const char * name) {
  xcint_assure_setup();
  int k;
  if ((k = xcint_lookup_functional(name)) >= 0)
    return xcint_funs[k].long_description;
  else if ((k = xcint_lookup_parameter(name)) >= 0)
    return xcint_params[k].description;
  else if ((k = xcint_lookup_alias(name)) >= 0)
    return xcint_aliases[k].description;
  else
    return 0;
}

XCFunctional::XCFunctional() {
  for (int i = 0; i < XC_NR_FUNCTIONALS; ++i)
    settings[i] = 0;
  for (int i = XC_NR_FUNCTIONALS; i < XC_NR_PARAMETERS_AND_FUNCTIONALS; ++i)
    settings[i] = xcint_params[i].default_value;
}

namespace xcfun {
XCFunctional * xcfun_new() {
  xcint_assure_setup();
  return new XCFunctional();
}

void xcfun_delete(XCFunctional * fun) {
  if (!fun)
    return;
  delete fun;
}

int xcfun_set(XCFunctional * fun, const char * name, double value) {
  xcint_assure_setup();
  int item;
  if ((item = xcint_lookup_functional(name)) >= 0) {
    fun->settings[item] += value;
    // Do not extend list if functional is active
    bool found = false;
    for (int i = 0; i < fun->nr_active_functionals; i++)
      if (fun->active_functionals[i] == &xcint_funs[item]) {
        found = true;
        break;
      }
    if (!found) {
      fun->active_functionals[fun->nr_active_functionals++] = &xcint_funs[item];
      fun->depends |= xcint_funs[item].depends;
    }
    return 0;
  } else if ((item = xcint_lookup_parameter(name)) >= 0) {
    fun->settings[item] = value;
    return 0;
  } else if ((item = xcint_lookup_alias(name)) >= 0) {
    for (int i = 0; i < MAX_ALIAS_TERMS; i++) {
      if (!xcint_aliases[item].terms[i].name)
        break;
      // FIXME: Do not weight parameters with value for aliases, but what about
      // EXX?
      if (xcfun_set(fun,
                    xcint_aliases[item].terms[i].name,
                    value * xcint_aliases[item].terms[i].weight) != 0) {
        fprintf(stderr, "Trying to set %s\n", xcint_aliases[item].terms[i].name);
        xcfun::die("Alias with unknown terms, fix aliases.cpp", item);
      }
    }
    return 0;
  }
  return -1;
}

int xcfun_get(const XCFunctional * fun, const char * name, double * value) {
  xcint_assure_setup();
  int item;
  if ((item = xcint_lookup_functional(name)) >= 0) {
    *value = fun->settings[item];
    return 0;
  } else if ((item = xcint_lookup_parameter(name)) >= 0) {
    *value = fun->settings[item];
    return 0;
  }
  return -1;
}

bool xcfun_is_gga(const XCFunctional * fun) { return (fun->depends & XC_GRADIENT); }

bool xcfun_is_metagga(const XCFunctional * fun) {
  return (fun->depends & (XC_LAPLACIAN | XC_KINETIC));
}

int xcfun_eval_setup(XCFunctional * fun,
                     xcfun_vars vars,
                     xcfun_mode mode,
                     int order) {
  // Check that vars are enough for the functional
  if ((fun->depends & xcint_vars[vars].provides) != fun->depends) {
    return xcfun::XC_EVARS;
  }
  if ((order < 0 || order > XCFUN_MAX_ORDER) ||
      (mode == XC_PARTIAL_DERIVATIVES && order > 4))
    return xcfun::XC_EORDER;
  if (mode == XC_POTENTIAL) {
    // GGA potential needs full laplacian
    if ((fun->depends & XC_GRADIENT) &&
        !(vars == XC_A_2ND_TAYLOR || vars == XC_A_B_2ND_TAYLOR ||
          vars == XC_N_2ND_TAYLOR || vars == XC_N_S_2ND_TAYLOR)) {
      return xcfun::XC_EVARS | xcfun::XC_EMODE;
    }
    // No potential for meta gga's
    if (fun->depends & (XC_LAPLACIAN | XC_KINETIC))
      return xcfun::XC_EMODE;
  }
  fun->mode = mode;
  fun->vars = vars;
  fun->order = order;
  return 0;
}

int xcfun_user_eval_setup(XCFunctional * fun,
                          const int order,
                          const unsigned int func_type,
                          const unsigned int dens_type,
                          const unsigned int mode_type,
                          const unsigned int laplacian,
                          const unsigned int kinetic,
                          const unsigned int current,
                          const unsigned int explicit_derivatives) {
  xcfun_vars vars = xcfun_which_vars(
      func_type, dens_type, laplacian, kinetic, current, explicit_derivatives);
  xcfun_mode mode = xcfun_which_mode(mode_type);
  return xcfun::xcfun_eval_setup(fun, vars, mode, order);
}

int xcfun_input_length(const XCFunctional * fun) {
  return xcint_vars[fun->vars].len;
}

int xcfun_output_length(const XCFunctional * fun) {
  if (fun->mode == XC_MODE_UNSET)
    xcfun::die("xc_output_length() called before a mode was succesfully set", 0);
  if (fun->vars == XC_VARS_UNSET)
    xcfun::die("xc_output_length() called before variables were succesfully set", 0);
  if (fun->order == -1)
    xcfun::die("xc_output_length() called before the order were succesfully set", 0);
  if (fun->mode == XC_PARTIAL_DERIVATIVES) {
    return taylorlen(xcint_vars[fun->vars].len, fun->order);
  } else if (fun->mode == XC_POTENTIAL) {
    if (fun->vars == XC_A || fun->vars == XC_A_2ND_TAYLOR)
      return 2; // Energy+potential
    else
      return 3; // Spin-resolved potential
  } else {
    xcfun::die("XC_CONTRACTED not implemented in xc_output_length()", 0);
    return 0;
  }
}

void xcfun_eval(const XCFunctional * fun, const double input[], double output[]) {
  if (fun->mode == XC_MODE_UNSET)
    xcfun::die("xc_eval() called before a mode was successfully set", 0);
  if (fun->vars == XC_VARS_UNSET)
    xcfun::die("xc_eval() called before variables were successfully set", 0);
  if (fun->order == -1 && fun->mode != XC_POTENTIAL)
    xcfun::die("xc_eval() called before the order was successfully set", 0);
  if (fun->mode == XC_PARTIAL_DERIVATIVES) {
    switch (fun->order) {
      case 0: {
        typedef ctaylor<ireal_t, 0> ttype;
        int inlen = xcint_vars[fun->vars].len;
        ttype in[XC_MAX_INVARS], out = 0;
        for (int i = 0; i < inlen; i++)
          in[i] = input[i];
        densvars<ttype> d(fun, in);
        for (int i = 0; i < fun->nr_active_functionals; i++)
          out += fun->settings[fun->active_functionals[i]->id] *
                 fun->active_functionals[i]->fp0(d);
        output[0] = out.get(CNST);
      } break;
#if XCFUN_MAX_ORDER >= 1
      case 1: {
        int inlen = xcint_vars[fun->vars].len;
        {
          typedef ctaylor<ireal_t, 2> ttype2;
          ttype2 in2[XC_MAX_INVARS], out2 = 0;
          for (int i = 0; i < inlen; i++)
            in2[i] = input[i];
          for (int j = 0; j < inlen / 2; j++) {
            in2[2 * j].set(VAR0, 1);
            in2[2 * j + 1].set(VAR1, 1);
            densvars<ttype2> d(fun, in2);
            out2 = 0;
            for (int i = 0; i < fun->nr_active_functionals; i++)
              out2 += fun->settings[fun->active_functionals[i]->id] *
                      fun->active_functionals[i]->fp2(d);
            in2[2 * j] = input[2 * j];
            in2[2 * j + 1] = input[2 * j + 1];
            output[2 * j + 1] = out2.get(VAR0); // First derivatives
            output[2 * j + 2] = out2.get(VAR1); // First derivatives
          }
          output[0] = out2.get(CNST); // Energy
        }
        if (inlen & 1) {
          typedef ctaylor<ireal_t, 1> ttype;
          int inlen = xcint_vars[fun->vars].len;
          ttype in[XC_MAX_INVARS], out = 0;
          for (int i = 0; i < inlen; i++)
            in[i] = input[i];
          int j = inlen - 1;
          {
            in[j].set(VAR0, 1);
            densvars<ttype> d(fun, in);
            out = 0;
            for (int i = 0; i < fun->nr_active_functionals; i++)
              out += fun->settings[fun->active_functionals[i]->id] *
                     fun->active_functionals[i]->fp1(d);
            in[j] = input[j];
            output[j + 1] = out.get(VAR0); // First derivatives
          }
          output[0] = out.get(CNST); // Energy
        }
      } break;
#endif
#if XCFUN_MAX_ORDER >= 2
#if XCFUN_MAX_ORDER >= 3
      // Do the third order derivatives here, then use the second order code. This is
      // getting expensive..
      case 3: {
        typedef ctaylor<ireal_t, 3> ttype;
        int inlen = xcint_vars[fun->vars].len;
        ttype in[XC_MAX_INVARS], out = 0;
        for (int i = 0; i < inlen; i++)
          in[i] = input[i];
        int k = 1 + inlen + (inlen * (inlen + 1)) / 2;
        for (int i = 0; i < inlen; i++) {
          in[i].set(VAR0, 1);
          for (int j = i; j < inlen; j++) {
            in[j].set(VAR1, 1);
            for (int s = j; s < inlen; s++) {
              in[s].set(VAR2, 1);
              densvars<ttype> d(fun, in);
              out = 0;
              for (int n = 0; n < fun->nr_active_functionals; n++)
                out += fun->settings[fun->active_functionals[n]->id] *
                       fun->active_functionals[n]->fp3(d);
              output[k++] = out.get(VAR0 | VAR1 | VAR2); // Third derivative
              in[s].set(VAR2, 0);
            }
            in[j].set(VAR1, 0);
          }
          in[i] = input[i];
        }
      }
#endif
      case 2: {
        typedef ctaylor<ireal_t, 2> ttype;
        int inlen = xcint_vars[fun->vars].len;
        ttype in[XC_MAX_INVARS], out = 0;
        for (int i = 0; i < inlen; i++)
          in[i] = input[i];
        int k = inlen + 1;
        for (int i = 0; i < inlen; i++) {
          in[i].set(VAR0, 1);
          for (int j = i; j < inlen; j++) {
            in[j].set(VAR1, 1);
            densvars<ttype> d(fun, in);
            out = 0;
            for (int n = 0; n < fun->nr_active_functionals; n++)
              out += fun->settings[fun->active_functionals[n]->id] *
                     fun->active_functionals[n]->fp2(d);
            output[k++] = out.get(VAR0 | VAR1); // Second derivative
            in[j].set(VAR1, 0);                 // slightly pessimized
          }
          output[i + 1] = out.get(VAR0); // First derivative
          in[i] = input[i];
        }
        output[0] = out.get(CNST); // Energy
      } break;
#endif
      default:
        xcfun::die("FIXME: Order too high for partial derivatives in xc_eval",
                   fun->order);
    }
  } else if (fun->mode == XC_CONTRACTED) {
#define DOEVAL(N, E)                                                                \
  if (fun->order == N) {                                                            \
    typedef ctaylor<ireal_t, N> ttype;                                              \
    int inlen = xcint_vars[fun->vars].len;                                          \
    ttype in[XC_MAX_INVARS], out = 0;                                               \
    int k = 0;                                                                      \
    for (int i = 0; i < inlen; i++)                                                 \
      for (int j = 0; j < (1 << fun->order); j++)                                   \
        in[i].set(j, input[k++]);                                                   \
    densvars<ttype> d(fun, in);                                                     \
    for (int i = 0; i < fun->nr_active_functionals; i++)                            \
      out += fun->settings[fun->active_functionals[i]->id] *                        \
             fun->active_functionals[i]->fp##N(d);                                  \
    for (int i = 0; i < (1 << fun->order); i++)                                     \
      output[i] = out.get(i);                                                       \
  } else
    FOR_EACH(XCFUN_MAX_ORDER, DOEVAL, )
    xcfun::die("bug! Order too high in XC_CONTRACTED", fun->order);
  } else if (fun->mode == XC_POTENTIAL) {
    // TODO: We shouldn't need the second density derivatives internally
    int inlen = xcint_vars[fun->vars].len;
    int npot; // One or two potentials
    int inpos = 0;
    if (inlen == 1 || inlen == 10)
      npot = 1;
    else {
      // nspin = 2 case. More complicated to find the
      // beta density as it depends on whether this is an lda or gga calculation.
      npot = 2;
      if (inlen == 2)
        inpos = 1;
      else if (inlen == 20)
        inpos = 10;
    }
    {
      typedef ctaylor<ireal_t, 1> ttype;
      ttype in[XC_MAX_INVARS], out = 0;
      for (int i = 0; i < inlen; i++)
        in[i] = input[i];
      for (int j = 0; j < npot; j++) {
        in[j * inpos].set(VAR0, 1);
        densvars<ttype> d(fun, in);
        out = 0;
        for (int i = 0; i < fun->nr_active_functionals; i++)
          out += fun->settings[fun->active_functionals[i]->id] *
                 fun->active_functionals[i]->fp1(d);
        in[j * inpos] = input[j * inpos];
        output[j + 1] = out.get(VAR0); // First derivatives
      }
      //            printf("Potentials 1 %f\n",output[1]);
      output[0] = out.get(CNST); // Energy
    }
    if (fun->depends & XC_GRADIENT) // GGA potential
    {
      /*
         v = dE/dn - nabla.dE/dg
       */
      typedef ctaylor<ireal_t, 2> ttype;
      ttype in[XC_MAX_INVARS];
      // n gx gy gz xx xy xz yy yz zz
      // 0 1  2  3  4  5  6  7  8  9
      if (fun->vars == XC_A_2ND_TAYLOR || fun->vars == XC_N_2ND_TAYLOR) {
        ttype out = 0;
        // d/dx
        in[0] = ttype(input[0], VAR0, input[1]);
        for (int i = 0; i < 3; i++)
          in[1 + i] = ttype(input[1 + i], VAR0, input[4 + i]);
        for (int i = 4; i < 10; i++)
          in[i] = 0;        // TODO: remove these vars, fun does not depend on them
        in[1].set(VAR1, 1); // d/dgx
        {
          densvars<ttype> d(fun, in);
          for (int i = 0; i < fun->nr_active_functionals; i++)
            out += fun->settings[fun->active_functionals[i]->id] *
                   fun->active_functionals[i]->fp2(d);
        }

        // d/dy
        in[0] = ttype(input[0], VAR0, input[2]);
        in[1] = ttype(input[1], VAR0, input[5]);
        in[2] = ttype(input[2], VAR0, input[7]);
        in[3] = ttype(input[3], VAR0, input[8]);
        in[2].set(VAR1, 1); // d/dgy
        {
          densvars<ttype> d(fun, in);
          for (int i = 0; i < fun->nr_active_functionals; i++)
            out += fun->settings[fun->active_functionals[i]->id] *
                   fun->active_functionals[i]->fp2(d);
        }
        // d/dz
        in[0] = ttype(input[0], VAR0, input[3]);
        in[1] = ttype(input[1], VAR0, input[6]);
        in[2] = ttype(input[2], VAR0, input[8]);
        in[3] = ttype(input[3], VAR0, input[9]);
        in[3].set(VAR1, 1); // d/dgz
        {
          densvars<ttype> d(fun, in);
          for (int i = 0; i < fun->nr_active_functionals; i++)
            out += fun->settings[fun->active_functionals[i]->id] *
                   fun->active_functionals[i]->fp2(d);
        }
        output[1] -= out.get(
            VAR0 | VAR1); // Subtract divergence of dE/dg from lda part of potential
      } else {
        // M Seth July-August 2011
        // Unrestricted GGA potential
        // a gx gy gz xx xy xz yy yz zz
        // 0 1  2  3  4  5  6  7  8  9
        // b  gx  gy  gz  xx  xy  xz  yy  yz  zz
        // 10 11  12  13  14  15  16  17  18  19
        {
          // j = 0 alpha and j = 1 beta
          for (int j = 0; j < 2; j++) {
            ttype out = 0;
            // Point to the correct set of values from the input
            int offset = 10;
            // d/dx
            in[0] = ttype(input[0], VAR0, input[1]);
            in[1] = ttype(input[1], VAR0, input[4]);
            in[2] = ttype(input[2], VAR0, input[5]);
            in[3] = ttype(input[3], VAR0, input[6]);
            in[0 + offset] = ttype(input[0 + offset], VAR0, input[1 + offset]);
            in[1 + offset] = ttype(input[1 + offset], VAR0, input[4 + offset]);
            in[2 + offset] = ttype(input[2 + offset], VAR0, input[5 + offset]);
            in[3 + offset] = ttype(input[3 + offset], VAR0, input[6 + offset]);
            for (int i = 4; i < 10; i++)
              in[i] = 0; // TODO: remove these vars, fun does not depend on them
            in[1 + offset * j].set(VAR1, 1); // d/dgx(a/b)
            {
              densvars<ttype> d(fun, in);
              for (int i = 0; i < fun->nr_active_functionals; i++)
                out += fun->settings[fun->active_functionals[i]->id] *
                       fun->active_functionals[i]->fp2(d);
            }
            // d/dy
            in[0] = ttype(input[0], VAR0, input[2]);
            in[1] = ttype(input[1], VAR0, input[5]);
            in[2] = ttype(input[2], VAR0, input[7]);
            in[3] = ttype(input[3], VAR0, input[8]);
            in[0 + offset] = ttype(input[0 + offset], VAR0, input[2 + offset]);
            in[1 + offset] = ttype(input[1 + offset], VAR0, input[5 + offset]);
            in[2 + offset] = ttype(input[2 + offset], VAR0, input[7 + offset]);
            in[3 + offset] = ttype(input[3 + offset], VAR0, input[8 + offset]);
            in[2 + offset * j].set(VAR1, 1); // d/dgy(a/b)
            {
              densvars<ttype> d(fun, in);
              for (int i = 0; i < fun->nr_active_functionals; i++)
                out += fun->settings[fun->active_functionals[i]->id] *
                       fun->active_functionals[i]->fp2(d);
            }
            // d/dz
            in[0] = ttype(input[0], VAR0, input[3]);
            in[1] = ttype(input[1], VAR0, input[6]);
            in[2] = ttype(input[2], VAR0, input[8]);
            in[3] = ttype(input[3], VAR0, input[9]);
            in[0 + offset] = ttype(input[0 + offset], VAR0, input[3 + offset]);
            in[1 + offset] = ttype(input[1 + offset], VAR0, input[6 + offset]);
            in[2 + offset] = ttype(input[2 + offset], VAR0, input[8 + offset]);
            in[3 + offset] = ttype(input[3 + offset], VAR0, input[9 + offset]);
            in[3 + offset * j].set(VAR1, 1); // d/dgy(a/b)
            {
              densvars<ttype> d(fun, in);
              for (int i = 0; i < fun->nr_active_functionals; i++)
                out += fun->settings[fun->active_functionals[i]->id] *
                       fun->active_functionals[i]->fp2(d);
            }
            output[j + 1] -= out.get(
                VAR0 |
                VAR1); // Subtract divergence of dE/dg from lda part of potential
          }
        }
        //	      xcfun::die("TODO: AB GGA potential",fun->mode);
      }
    }
  } else {
    xcfun::die("Illegal mode in xc_eval()", fun->mode);
  }
}

void xcfun_eval_vec(const XCFunctional * fun,
                    int nr_points,
                    const double density[],
                    int density_pitch,
                    double result[],
                    int result_pitch) {
  for (int i = 0; i < nr_points; i++)
    xcfun_eval(fun, density + i * density_pitch, result + i * result_pitch);
}
} // namespace xcfun

xcfun_t * xcfun_new() { return AS_TYPE(xcfun_t, xcfun::xcfun_new()); }

void xcfun_delete(xcfun_t * fun) { xcfun::xcfun_delete(AS_TYPE(XCFunctional, fun)); }

int xcfun_set(xcfun_t * fun, const char * name, double value) {
  return xcfun::xcfun_set(AS_TYPE(XCFunctional, fun), name, value);
}

int xcfun_get(const xcfun_t * fun, const char * name, double * value) {
  return xcfun::xcfun_get(AS_CTYPE(XCFunctional, fun), name, value);
}

bool xcfun_is_gga(const xcfun_t * fun) {
  return xcfun::xcfun_is_gga(AS_CTYPE(XCFunctional, fun));
}

bool xcfun_is_metagga(const xcfun_t * fun) {
  return xcfun::xcfun_is_metagga(AS_CTYPE(XCFunctional, fun));
}

int xcfun_eval_setup(xcfun_t * fun, xcfun_vars vars, xcfun_mode mode, int order) {
  return xcfun::xcfun_eval_setup(AS_TYPE(XCFunctional, fun), vars, mode, order);
}

int xcfun_user_eval_setup(xcfun_t * fun,
                          const int order,
                          const unsigned int func_type,
                          const unsigned int dens_type,
                          const unsigned int mode_type,
                          const unsigned int laplacian,
                          const unsigned int kinetic,
                          const unsigned int current,
                          const unsigned int explicit_derivatives) {
  return xcfun::xcfun_user_eval_setup(AS_TYPE(XCFunctional, fun),
                                      order,
                                      func_type,
                                      dens_type,
                                      mode_type,
                                      laplacian,
                                      kinetic,
                                      current,
                                      explicit_derivatives);
}

int xcfun_input_length(const xcfun_t * fun) {
  return xcfun::xcfun_input_length(AS_CTYPE(XCFunctional, fun));
}

int xcfun_output_length(const xcfun_t * fun) {
  return xcfun::xcfun_output_length(AS_CTYPE(XCFunctional, fun));
}

void xcfun_eval(const xcfun_t * fun, const double input[], double output[]) {
  return xcfun::xcfun_eval(AS_CTYPE(XCFunctional, fun), input, output);
}

void xcfun_eval_vec(const xcfun_t * fun,
                    int nr_points,
                    const double density[],
                    int density_pitch,
                    double result[],
                    int result_pitch) {
  xcfun::xcfun_eval_vec(AS_CTYPE(XCFunctional, fun),
                        nr_points,
                        density,
                        density_pitch,
                        result,
                        result_pitch);
}
