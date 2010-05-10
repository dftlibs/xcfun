#include "functional.h"
#include "xcfun_internal.h"
#include <cstring>
#include <cstdio>


functional::functional()
{
  for (int i=0;i<=XC_MAX_NVAR;i++)
    for (int j=0;j<=XC_MAX_ORDER;j++)
      ftab[i][j] = 0;
  m_name = 0;
  m_oneliner = 0;
  m_description = 0;
  test_mode = -1;
}

functional::~functional()
{
  if (test_mode != -1)
    {
      delete[] test_input;
      delete[] test_output;
    }
}



void functional::describe(const char *name, int type, const char *oneliner,
			  const char *reference)
{
  assert(type >= 0 and type < XC_NR_TYPES);
  m_type = type;
  m_name = name;
  m_oneliner = oneliner;
  m_description = reference;
}

void functional::parameter(const char *name, const char *description, 
			   double default_value)
{
  xc_get_settings().new_setting(name,default_value,description,description);
}

void functional::add_test(int mode, int order,
		      const double *test_in,
		      const double *test_out,
		      double threshold)
{
  assert(test_mode == -1 && "Multiple tests not yet implemeted");
  test_mode = mode;
  test_order = order;

  test_input = new double[xc_input_length(mode,m_type)];
  for (int i=0;i<xc_input_length(mode,m_type);i++)
    test_input[i] = test_in[i];

  test_output = new double[xc_output_length(mode,m_type,order)];
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
  array<functional *> a = xc_get_functional_array();
  for (int i=0;i<a.size();i++)
    if (strcmp(a[i]->m_name,name) == 0)
      return a[i];
  return 0;
}

// Create a new functional object in the internal database,
// return a reference to this object (to be filled with data)
functional &xc_new_functional(void)
{
  functional *f = new functional;
  xc_get_functional_array().push_back(f);
  return *f;
}

int xc_run_functional_setup(void (*setup)(functional &))
{
  functional &f = xc_new_functional();
  setup(f);
  xc_get_settings().new_setting(f.m_name,0,f.m_oneliner,f.m_description);
  return f.validate();
}

