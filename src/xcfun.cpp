
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "xcfun_internal.h"

static bool xcfun_is_setup = false;

static void xcint_assure_setup()
{
  if (!xcfun_is_setup)
    {
      xcint_setup_functionals();
      xcfun_is_setup = true;
    }
}

void xc_die(const char *message, int code)
{
  fprintf(stderr,"XCFun fatal error %i: ",code);
  fprintf(stderr,"%s",message);
  fprintf(stderr,"\n");
  exit(-1);
}

extern "C"
double xcfun_version(void)
{
  return 0.9;
}

extern "C"
const char *xcfun_splash(void)
{
  return "XCFun DFT library, Copyright 2009-2010 Ulf Ekstrom\n"
    "and contributors. See http://admol.org/xcfun for more information.\n"
    "This is free software; see the source code for copying conditions.\n"
    "There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or\n"
    "FITNESS FOR A PARTICULAR PURPOSE. For details see the documentation.\n";
}

int xc_input_length(xc_functional fun, int mode, int type)
{
  static int tab[XC_NR_MODES][XC_NR_TYPES] = 
    {{1,2,3},{1,2,3},{2,5,7},{2,5,7}};
  assert(mode>=0 and mode <= XC_NR_MODES);
  assert(type>=0 and type <= XC_NR_TYPES);
  return tab[mode][type];
}

int xc_output_length(int mode, int type, int order)
{
  return taylorlen(xc_input_length(mode,type),order);
}

settings_database &xc_get_settings()
{
  static settings_database b;
  return b;
}

setting xc_param_lookup(const xc_functional::xc_functional_data *params, 
			const char *name)
{
  return params->settings->lookup(name);
}

double xc_param_get(const xc_functional::xc_functional_data *params,
		    const setting &s)
{
  return params->settings->get(s);
}


xc_functional_data::initialize(void)
{
  xcint_assure_setup();
  mode = XC_VARS_AB; //Use alpha/beta as default
  type = XC_LDA; // LDA is smallest type, so assume this initially
  max_order = -1;
  // Just default settings to start with
  settings = xc_get_settings().new_user_settings();
}

xc_functional_data::destroy(void)
{
  delete settings;
}

void xc_functional_data::regularize_density(double *density)
{
  //TODO: Actually do something here.
}

int xc_functional_data::get_type(void) const
{
  return type;
}
void xc_functional_data::set_mode(int mode)
{
  if (!(mode>=0 and mode < XC_NR_MODES))
    xc_die("Invalid mode to xc_functional::set_mode()",mode);
  this->mode = mode;
  find_max_order();
}

void xc_functional_data::find_max_order(void)
{
  max_order = -1;
  while (max_order < XC_MAX_ORDER and xc_evaluator_lookup(mode,type,max_order+1))
    max_order++;
}

int xc_functional_data::get_max_order(void) const
{
  return max_order;
}


int xc_functional_data::input_length(void) const
{
  return xc_input_length(mode,type);
}

int xc_functional_data::output_length(int order) const
{
  return taylorlen(input_length(),order);
}

int xc_functional_data::
derivative_index(const int exponents[]) const
{
  int nvar = input_length();
  int N = 0;
  for (int i=0;i<nvar;i++)
    N += exponents[i];
  int i = 0, idx = 0;
  N--;
  while (N >= 0)
    {
      idx += taylorlen(nvar-i,N);
      N -= exponents[i];
      i++;
    }
  return idx;
}


// API starts here
extern "C"
xc_functional xc_new_functional(void)
{
  xc_functional_data *p = malloc(sizeof*p);
  p->initialize();
  return p;
}

extern "C"
void xc_free_functional(xc_functional fun)
{
  fun->destroy();
  free(fun);
}

void xc_regularize_density(xc_functional fun,double *density)
{
  fun->regularize_density(density);
}

void xc_eval(xc_functional fun, int order, int nr_points,
	     const double *density, double *result)
{
  evaluator ev = xc_evaluator_lookup(fun->mode,fun->type,order);
  if (!ev)
    {
      fprintf(stderr,"XCFun error in eval()\n");
      fprintf(stderr,"mode: %i\n",fun->mode);
      fprintf(stderr,"type: %i\n",fun->type);
      xc_die("eval(): Functional not available for order",order);
    }
  int inlen = xc_input_length(fun);
  int outlen = xc_output_length(fun,order);
  for (int i=0;i<nr_points;i++)
    ev(*fun,result+i*outlen,density+i*inlen);
}

int xc_get_type(xc_functional fun)
{
  return fun->get_type();
}

int xc_max_order(xc_functional fun)
{
  return fun->get_max_order();
}

int xc_input_length(xc_functional fun)
{
  return fun->input_length();
}

int xc_output_length(xc_functional fun, int order)
{
  return fun->output_length(order);
}

int xc_derivative_index(xc_functional fun, const int derivative[])
{
  return fun->derivative_index(derivative);
}

void xc_set_mode(xc_functional fun, int mode)
{
  fun->set_mode(mode);
}

int xc_set_setting(xc_functional fun, int setting_nr, double value)
{
  if (is_functional(name))
    {
      functional *fun = xc_get_functional_by_name(name);
      assert(fun);
      if (d->type < fun->m_type)
	d->type = fun->m_type;
      int i = d->active_functionals.find(fun);
      if (i>=0)
	{
	  d->weights[i] = value;
	}
      else
	{
	  d->active_functionals.push_back(fun);
	  d->weights.push_back(value);
	}
    }
  return d->settings->set(name,value);
}

double get_setting(const char *name) const
{
  return d->settings->get(name);
}

bool is_functional(const char *name) const
{
  array<functional *> &a = xc_get_functional_array();
  for (int i=0;i<a.size();i++)
    if (strcmp(name,a[i]->m_name) == 0)
      return true;
  return false;
}

bool is_set(const char *name) const
{
  return d->settings->is_set(name);
}

const char *setting_name(int n) const
{
  if ( n<0 or n >= d->settings->nr_settings() )
    return 0;
  else
    return d->settings->setting_name(n);
}

const char *setting_short_description(const char *name) const
{
  int i = d->settings->index_of(name);
  if (i >= 0)
    return d->settings->setting_short_description(i);
  else
    return 0;
}

const char *setting_long_description(const char *name) const
{
  int i = d->settings->index_of(name);
  if (i >= 0)
    return d->settings->setting_long_description(i);
  else
    return 0;
}

static int run_tests(functional *fun)
{
  if (fun->test_mode == -1)
    return -1;
  int n = xc_output_length(fun->test_mode,fun->m_type,fun->test_order);
  double *out = new double[n];
  double *reference = fun->test_output;
  xc_functional xf;
  xf.set_mode(fun->test_mode);
  xf.set_setting(fun->m_name,1.0); //Weight 1.0
  xf.eval(out,fun->test_order,fun->test_input);
  int nerr = 0;
  for (int i=0;i<n;i++)
    if (fabs(out[i] - fun->test_output[i]) > 
	fabs(fun->test_output[i]*fun->test_threshold))
      nerr++;

  if (nerr > 0)
    {
      fprintf(stderr,"Error detected in functional %s with tolerance %g:\n",
	      fun->m_name,fun->test_threshold);
      fprintf(stderr,"Abs.Error \tComputed              Reference\n");
      for (int i=0;i<n;i++)
	{
	  fprintf(stderr,"%.1e",fabs(out[i]-reference[i]));
	  fprintf(stderr,"    %+.16f \t%+.16f",out[i],reference[i]);
	  if (fabs(out[i] - reference[i]) > 
	      fabs(reference[i]*fun->test_threshold))
	     fprintf(stderr," *");
	  fprintf(stderr,"\n");
	}
    }
  return nerr;
}

int xcfun_test(void)
{
  int nr_failed = 0;
  array<functional *> &a = xc_get_functional_array();
  for (int i=0;i<a.size();i++)
    {
      int res = run_tests(a[i]);
      if (res != 0 and res != -1)
	nr_failed++;
    }
  return nr_failed;
}
