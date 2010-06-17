#include "functional.h"
#include "xcfun_internal.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

void functional::construct()
{
  for (int i=0;i<=XC_MAX_NVAR;i++)
    for (int j=0;j<=XC_MAX_ORDER;j++)
      ftab[i][j] = 0;
  m_name = XC_NR_PARAMS; // Must be set to a smaller value before use!
  test_mode = -1;
}

void functional::destroy()
{
  if (test_mode != -1)
    {
      free(test_input);
      free(test_output);
    }
}

void functional::describe(enum xc_parameters weight_param, int type, 
			  const char *oneliner,
			  const char *reference)
{
  assert(type >= 0 and type < XC_NR_TYPES);
  m_type = type;
  m_name = weight_param;
  xcint_set_short_description(weight_param, oneliner);
  xcint_set_long_description(weight_param, reference);
}

void functional::describe_parameter(enum xc_parameters param, const char *description, 
			   double default_value)
{
  xcint_set_short_description(param, description);
  xcint_set_default(param,default_value);
}

void functional::add_test(int mode, int order,
		      const double *test_in,
		      const double *test_out,
		      double threshold)
{
  assert(test_mode == -1 && "Multiple tests not yet implemeted");
  test_mode = mode;
  test_order = order;

  test_input = (double *)malloc(xcint_input_length(mode,m_type)*sizeof*test_input);
  for (int i=0;i<xcint_input_length(mode,m_type);i++)
    test_input[i] = test_in[i];

  test_output = (double *)malloc(xcint_output_length(mode,m_type,order)*sizeof*test_output);
  for (int i=0;i<xcint_output_length(mode,m_type,order);i++)
    test_output[i] = test_out[i];

  test_threshold = threshold;
}

int functional::validate()
{
  //TODO: Check that the oneliner has no \n
  if (test_mode == -1)
    {
      fprintf(stderr,"XCFUN WARNING: Functional %s has no test\n",xc_name(m_name));
      return -1;
    }
  for (int i=0;i<=XC_MAX_NVAR;i++)
    for (int j=0;j<=XC_MAX_ORDER;j++)
      if (ftab[i][j])
	return 0;
  fprintf(stderr,"XCFUN WARNING: Functional %s has no implementation\n",xc_name(m_name));
  return -1;
}

functional *xc_functional_by_name(const char *name)
{
  for (int i=0;i<XC_NR_PARAMS;i++)
    if (xcint_functional(i) and
	strcmp(xc_name(i),name) == 0)
      return xcint_functional(i);
  return 0;
}

int xc_run_functional_setup(void (*setup)(functional &))
{
  functional *f = (functional *)malloc(sizeof(functional));
  f->construct();
  setup(*f);
  xcint_set_functional(f->m_name,f);
  return f->validate();
}

