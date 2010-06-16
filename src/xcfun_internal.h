#ifndef XCFUN_INTERNAL_H
#define XCFUN_INTERNAL_H

#include "xcfun.h"
#include "functional.h"
#include "array.h"
#include "config.h"
#include "settings.h"

struct xc_functional_data
{
public:
  void initialize(void);
  void destroy(void);

  void set_mode(int mode);

  int get_type(void) const;
  int get_max_order(void) const;

  void regularize_density(double *density);
  void eval(double *res, int order, const double *dens) const;

  bool is_functional(const char *name);

  int input_length(void) const;
  int output_length(int order) const;
  int derivative_index(const int derivative[]) const;

  void find_max_order(void);

  int mode; // One of XC_VARS_*
  int type; // LDA, GGA etc
  int max_order; // Maximum derivative order with current settings

  user_settings *settings;
  array<functional *> active_functionals;
  array<double> weights;
};

void xc_die(const char *message, int code);
int xc_input_length(int mode, int type);
int xc_output_length(int mode, int type, int order);

settings_database &xc_get_settings();
typedef void (*evaluator)(const xc_functional::xc_functional_data &fun, double *, const double *);
evaluator xc_evaluator_lookup(int mode, int type, int order);

void xcint_setup_functionals();

#endif
