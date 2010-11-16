#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

#include "xcint.h"
#include "taylor.h"
#include "config.h"
#include "specmath.h"

typedef ireal_t parameter;

// Variables for expressing functionals, these are redundant because
// different functionals have different needs.
template<class T>
struct densvars
{
  const xc_functional_obj *parent;
  double get_param(enum xc_parameter p) const
  {
    return parent->settings[p];
  }

  
  T a, b, gaa, gab, gbb;
  /* na+nb, na-nb, (grad n)^2, (grad n).(grad s), (grad s)^2 */
  T n, s, gnn, gns, gss; 

  T tau, taua, taub; // Kinetic energy densities.

  T lapa,lapb; // Density Laplacians

  T zeta; //s/n
  T r_s; // (3/4pi)^1/3*n^(-1/3)
  T n_m13; // pow(n,-1.0/3.0)
  T a_43, b_43; // pow(a,4.0/3.0), pow(b,4.0/3.0)
};

template<int F>
struct functional
{
  static const char *name;
  static const char *short_description;
  static const char *long_description;
  static const double test_in[];
  static const double test_out[];
  static const double test_threshold;
  static const enum xc_vars test_vars;
  static const enum xc_mode test_mode;
  static const int test_order;
};

template<int F, class T>
struct erg
{
  typedef T (*t)(const densvars<T> &d);
  static T (*f)(const densvars<T> &d);
};

// By default we don't have any implementations for a functional
// template<int F, class T>  T (*erg<F,T>::f)(const densvars<T> &) = 0;

#define ENERGY_FUNCTION(FUN) ef_##FUN

#define IN(FUN,NVAR,NDEG) template<> erg<FUN,taylor<ireal_t,NVAR,NDEG> >::t erg<FUN, taylor<ireal_t,NVAR,NDEG> >::f = ef_##FUN;
#define DU(FUN,NVAR,NDEG) template<> erg<FUN,taylor<ireal_t,NVAR,NDEG> >::t erg<FUN, taylor<ireal_t,NVAR,NDEG> >::f = 0;

#define ON(FUN,NVAR) IN(FUN,NVAR,0) IN(FUN,NVAR,1) IN(FUN,NVAR,2)
#define OFF(FUN,NVAR) DU(FUN,NVAR,0) DU(FUN,NVAR,1) DU(FUN,NVAR,2)

#define NEW_LDA_FUNCTIONAL(FUN) template<> const char *functional<FUN>::name = #FUN;\
  ON(FUN,1) ON(FUN,2) OFF(FUN,3) OFF(FUN,5) OFF(FUN,7)

#define NEW_GGA_FUNCTIONAL(FUN) template<> const char *functional<FUN>::name = #FUN;\
  OFF(FUN,1) ON(FUN,2) OFF(FUN,3) ON(FUN,5) OFF(FUN,7)

#define NEW_TMGGA_FUNCTIONAL(FUN) template<> const char *functional<FUN>::name = #FUN;\
  OFF(FUN,1) OFF(FUN,2) ON(FUN,3) OFF(FUN,5) ON(FUN,7)

#define NEW_LTMGGA_FUNCTIONAL(FUN) template<> const char *functional<FUN>::name = #FUN;\
  OFF(FUN,1) OFF(FUN,2) OFF(FUN,3) OFF(FUN,5) OFF(FUN,7)

#define SHORT_DESCRIPTION(FUN) template<> const char *functional<FUN>::short_description
#define LONG_DESCRIPTION(FUN) template<> const char *functional<FUN>::long_description
#define TEST_IN(FUN)    template<> const double functional<FUN>::test_in[]
#define TEST_OUT(FUN)   template<> const double functional<FUN>::test_out[]
#define TEST_VARS(FUN)  template<> const enum xc_vars functional<FUN>::test_vars
#define TEST_MODE(FUN)  template<> const enum xc_mode functional<FUN>::test_mode
#define TEST_ORDER(FUN) template<> const int functional<FUN>::test_order
#define TEST_THRESHOLD(FUN) template<> const double functional<FUN>::test_threshold
#define NO_TEST(FUN) template<> TEST_IN(FUN) = {0}; TEST_OUT(FUN) = {0}; TEST_VARS(FUN) = XC_VARS_UNSET; TEST_MODE(FUN) = XC_MODE_UNSET; TEST_ORDER(FUN) = -1; TEST_THRESHOLD(FUN) = 0

#endif
