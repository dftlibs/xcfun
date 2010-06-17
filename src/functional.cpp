#include "functional.h"
#include "xcfun_internal.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

array<functional *> all_functionals;

functional::construct()
{
  for (int i=0;i<=XC_MAX_NVAR;i++)
    for (int j=0;j<=XC_MAX_ORDER;j++)
      ftab[i][j] = 0;
  m_name = 0;
  m_oneliner = 0;
  m_description = 0;
  test_mode = -1;
}

functional::destroy()
{
  if (test_mode != -1)
    {
      free(test_input);
      free(test_output);
    }
}



void functional::describe(enum xc_parameter weight_param, int type, 
			  const char *oneliner,
			  const char *reference)
{
  assert(type >= 0 and type < XC_NR_TYPES);
  m_type = type;
  m_name = weight_param;
  m_oneliner = oneliner;
  m_description = reference;
  xc_param_set_short_description(weight_param, oneliner);
  xc_param_set_long_description(weight_param, reference);
}

void functional::describe_parameter(enum xc_parameter param, const char *description, 
			   double default_value)
{
  xc_param_set_short_description(param, description);
}

void functional::add_test(int mode, int order,
		      const double *test_in,
		      const double *test_out,
		      double threshold)
{
  assert(test_mode == -1 && "Multiple tests not yet implemeted");
  test_mode = mode;
  test_order = order;

  test_input = malloc(xc_input_length(mode,m_type)*sizeof*test_input);
  for (int i=0;i<xc_input_length(mode,m_type);i++)
    test_input[i] = test_in[i];

  test_output = malloc(xc_output_length(mode,m_type,order)*sizeof*test_output);
  for (int i=0;i<xc_output_length(mode,m_type,order);i++)
    test_output[i] = test_out[i];

  test_threshold = threshold;
}

int functional::validate()
{
  if (!m_name)
    {
      fprintf(stderr,"XCFUN WARNING: Functional not given a name\n");
      return -1;
    }
  //TODO: Check that the oneliner has no \n
  if (test_mode == -1)
    {
      fprintf(stderr,"XCFUN WARNING: Functional %s has no test\n",m_name);
      return -1;
    }
  for (int i=0;i<=XC_MAX_NVAR;i++)
    for (int j=0;j<=XC_MAX_ORDER;j++)
      if (ftab[i][j])
	return 0;
  fprintf(stderr,"XCFUN WARNING: Functional %s has no implementation\n",m_name);
  return -1;
}


array<functional *> &xc_get_functional_array(void)
{
  static array<functional *> info;
  return info;
}

functional *xc_get_functional_by_name(const char *name)
{
  for (int i=0;i<all_functionals.size();i++)
    if (strcmp(all_functionals[i]->m_name,name) == 0)
      return all_functionals[i];
  return 0;
}

// Create a new functional object in the internal database,
// return a reference to this object (to be filled with data)
functional &xc_new_functional(void)
{
  functional *f = malloc(sizeof(functional));
  f->construct();
  all_functionals.push_back(f);
  return *f;
}

int xc_run_functional_setup(void (*setup)(functional &))
{
  functional &f = xc_new_functional();
  setup(f);
  return f.validate();
}

