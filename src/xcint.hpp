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

#include <cstdio>

#include "config.hpp"
#include "ctaylor.hpp"
#include "functionals/list_of_functionals.hpp"
#include "taylor.hpp"
#include "xcfun.h"

#define XC_MAX_ALIASES 60
#define MAX_ALIAS_TERMS 10
#define XC_MAX_INVARS 20

// Macros to iterate up to XC_MAX_ORDER
#define REP0(F, E) F(0, E)
#define REP1(F, E) REP0(F, E) F(1, E)
#define REP2(F, E) REP1(F, E) F(2, E)
#define REP3(F, E) REP2(F, E) F(3, E)
#define REP4(F, E) REP3(F, E) F(4, E)
#define REP5(F, E) REP4(F, E) F(5, E)
#define REP6(F, E) REP5(F, E) F(6, E)
#define REP7(F, E) REP6(F, E) F(7, E)
#define REP8(F, E) REP7(F, E) F(8, E)
#define REP9(F, E) REP8(F, E) F(9, E)
#define XFOR_EACH(N, F, E) REP##N(F, E)
#define FOR_EACH(N, F, E) XFOR_EACH(N, F, E)

#define XC_DENSITY 1
#define XC_GRADIENT 2
#define XC_LAPLACIAN 4
#define XC_KINETIC 8
#define XC_JP 16

typedef void (*evaluator)(struct xc_functional_obj *, const double *, double *);

struct functional_data;

struct xc_functional_obj {
  double settings[XC_NR_PARAMETERS_AND_FUNCTIONALS];
  int nr_active_functionals;
  const functional_data * active_functionals[XC_NR_FUNCTIONALS];
  enum xc_mode mode;
  enum xc_vars vars;
  int order;
  int depends; // XC_DENSITY, gradient etc
};

template <typename T> struct densvars;

struct functional_data {
  const char * short_description;
  const char * long_description;
  int depends; // XC_DENSITY | XC_GRADIENT etc
#define FP(N, E) ctaylor<ireal_t, N> (*fp##N)(const densvars<ctaylor<ireal_t, N>> &);
  FOR_EACH(XC_MAX_ORDER, FP, )
  enum xc_vars test_vars;
  enum xc_mode test_mode;
  int test_order;
  double test_threshold;
  double test_in[16]; // Increase dimensions if future tests require it
  double test_out[128];

  enum xc_functional_id id;
  const char * name; // Set up automatically from the symbol
};

typedef ireal_t parameter;

struct parameter_data {
  const char * description;
  parameter default_value;

  const char * name; // Set up automatically
};

struct vars_data {
  const char * symbol;
  int len;
  int provides; // XC_DENSITY | XC_GRADIENT etc
};

struct alias_data {
  const char * name;
  const char * description;
  struct {
    const char * name;
    double weight;
  } terms[MAX_ALIAS_TERMS];
};

extern functional_data xcint_funs[XC_NR_FUNCTIONALS];
extern parameter_data xcint_params[XC_NR_PARAMETERS_AND_FUNCTIONALS];
extern vars_data xcint_vars[XC_NR_VARS];
extern alias_data * xcint_aliases;

void xcint_assure_setup();

evaluator xcint_get_evaluator(enum xc_mode mode, enum xc_vars vars, int order);
extern "C" void xcint_die(const char * message, int code);
int xcint_write_fortran_module();

#if 0
static inline int taylorlen(int nvar, int ndeg)
{
  int len = 1;
  for (int k=1;k<=nvar;k++)
    {
      len *= ndeg + k;
      len /= k;
    }
  return len;
}
#endif

// The lookup functions return -1 if not found
// TODO: Case insensitive string comparison should be used
int xcint_lookup_functional(const char * name);
int xcint_lookup_parameter(const char * name);
int xcint_lookup_alias(const char * name);

#include "densvars.hpp"

// This gets filled in by the functional implementations
template <int FUN> struct fundat_db {
  static const char * symbol;
  static functional_data d;
};

template <int FUN> struct pardat_db {
  static const char * symbol;
  static parameter_data d;
};
