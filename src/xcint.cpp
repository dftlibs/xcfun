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

#include "xcint.hpp"
#include <cstdio>
#include <cstdlib>
#include <string.h>

#ifndef FORTRAN_INT
typedef int fortran_int_t;
#else
typedef FORTRAN_INT fortran_int_t;
#endif

functional_data xcint_funs[XC_NR_FUNCTIONALS];
parameter_data xcint_params[XC_NR_PARAMETERS_AND_FUNCTIONALS];

int xcint_lookup_functional(const char * name) {
  for (int i = 0; i < XC_NR_FUNCTIONALS; i++)
    if (strcasecmp(name, xcint_funs[i].name) == 0)
      return i;
  return -1;
}

int xcint_lookup_parameter(const char * name) {
  for (int i = XC_NR_FUNCTIONALS; i < XC_NR_PARAMETERS_AND_FUNCTIONALS; i++)
    if (strcasecmp(name, xcint_params[i].name) == 0)
      return i;
  return -1;
}

int xcint_lookup_alias(const char * name) {
  for (int i = 0; i < XC_MAX_ALIASES and xcint_aliases[i].name; i++)
    if (strcasecmp(name, xcint_aliases[i].name) == 0)
      return i;
  return -1;
}

template <int FUN> void xcint_functional_setup_helper() {
  if (!(fundat_db<FUN>::symbol[0] == 'X' and fundat_db<FUN>::symbol[1] == 'C' and
        fundat_db<FUN>::symbol[2] == '_'))
    xcint_die("Functional symbol does not start with XC_", FUN);
  fundat_db<FUN>::d.name = fundat_db<FUN>::symbol + 3;
  fundat_db<FUN>::d.id = (enum xc_functional_id)FUN;
  xcint_funs[FUN] = fundat_db<FUN>::d;
  xcint_functional_setup_helper<FUN + 1>();
}

template <> void xcint_functional_setup_helper<XC_NR_FUNCTIONALS>() {}

template <int FUN, int SPAN> struct retarded_helper {
  static void doit() {
    retarded_helper<FUN, SPAN / 2>::doit();
    retarded_helper<FUN + SPAN / 2, (SPAN + 1) / 2>::doit();
  }
};

template <int FUN> struct retarded_helper<FUN, 1> {
  static void doit() {
    if (!(fundat_db<FUN>::symbol[0] == 'X' and fundat_db<FUN>::symbol[1] == 'C' and
          fundat_db<FUN>::symbol[2] == '_'))
      xcint_die("Functional symbol does not start with XC_", FUN);
    fundat_db<FUN>::d.name = fundat_db<FUN>::symbol + 3;
    fundat_db<FUN>::d.id = (enum xc_functional_id)FUN;
    xcint_funs[FUN] = fundat_db<FUN>::d;
  }
};

template <int P> void xcint_parameter_setup_helper() {
  if (!(pardat_db<P>::symbol[0] == 'X' and pardat_db<P>::symbol[1] == 'C' and
        pardat_db<P>::symbol[2] == '_'))
    xcint_die("Symbol does not start with XC_", P);
  pardat_db<P>::d.name = pardat_db<P>::symbol + 3;
  xcint_params[P] = pardat_db<P>::d;
  xcint_parameter_setup_helper<P + 1>();
}

template <> void xcint_parameter_setup_helper<XC_NR_PARAMETERS_AND_FUNCTIONALS>() {}

vars_data xcint_vars[XC_NR_VARS] = {
    {"XC_A", 1, XC_DENSITY},
    {"XC_N", 1, XC_DENSITY},
    {"XC_A_B", 2, XC_DENSITY},
    {"XC_N_S", 2, XC_DENSITY},
    {"XC_A_GAA", 2, XC_DENSITY | XC_GRADIENT},
    {"XC_N_GNN", 2, XC_DENSITY | XC_GRADIENT},
    {"XC_A_B_GAA_GAB_GBB", 5, XC_DENSITY | XC_GRADIENT},
    {"XC_N_S_GNN_GNS_GSS", 5, XC_DENSITY | XC_GRADIENT},
    {"XC_A_GAA_LAPA", 3, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_A_GAA_TAUA", 3, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_N_GNN_LAPN", 3, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_GNN_TAUN", 3, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_A_B_GAA_GAB_GBB_LAPA_LAPB", 7, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_A_B_GAA_GAB_GBB_TAUA_TAUB", 7, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_N_S_GNN_GNS_GSS_LAPN_LAPS", 7, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_S_GNN_GNS_GSS_TAUN_TAUS", 7, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB",
     9,
     XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN | XC_KINETIC},
    {"XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB",
     11,
     XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN | XC_KINETIC | XC_JP},
    {"XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS",
     9,
     XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN | XC_KINETIC},
    {"XC_A_AX_AY_AZ", 4, XC_DENSITY | XC_GRADIENT},
    {"XC_A_B_AX_AY_AZ_BX_BY_BZ", 8, XC_DENSITY | XC_GRADIENT},
    {"XC_N_NX_NY_NZ", 4, XC_DENSITY | XC_GRADIENT},
    {"XC_N_S_NX_NY_NZ_SX_SY_SZ", 8, XC_DENSITY | XC_GRADIENT},
    {"XC_A_AX_AY_AZ_TAUA", 5, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB",
     10,
     XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_N_NX_NY_NZ_TAUN", 5, XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS",
     10,
     XC_DENSITY | XC_GRADIENT | XC_KINETIC},
    {"XC_A_2ND_TAYLOR", 10, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_A_B_2ND_TAYLOR", 20, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_2ND_TAYLOR", 10, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
    {"XC_N_S_2ND_TAYLOR", 20, XC_DENSITY | XC_GRADIENT | XC_LAPLACIAN},
};

void xcint_assure_setup() {
  static bool is_setup = false;
  if (!is_setup) {
    retarded_helper<0, XC_NR_FUNCTIONALS>::doit();
    //      xcint_functional_setup_helper<0>();
    xcint_parameter_setup_helper<XC_NR_FUNCTIONALS>();
#ifndef NDEBUG
    /* Verify that the variable definition is consistent. */
    for (int i = 0; i < XC_NR_VARS; i++) {
      assert(xcint_vars[i].len <= XC_MAX_INVARS);
    }
#endif
    is_setup = true;
  }
}

void xcint_die(const char * message, int code) {
  fprintf(stderr, "XCFun fatal error %i: ", code);
  fprintf(stderr, "%s", message);
  fprintf(stderr, "\n");
  exit(-1);
}
