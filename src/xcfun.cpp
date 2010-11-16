#include <cstdio>
#include <cstdlib>
#include "xcint.h"


xc_functional xc_new_functional()
{
  xc_functional fun = (xc_functional)malloc(sizeof*fun);
  if (!fun)
    xcint_die("Out of memory in xc_new_functional",0);
  fun->mode = XC_MODE_UNSET;
  fun->vars = XC_VARS_UNSET;
  fun->order = -1;
  for (int i=0;i<XC_NR_FUNCTIONALS;i++)
    fun->settings[i] = 0;
  for (int i=XC_NR_FUNCTIONALS;i<XC_NR_PARAMETERS_AND_FUNCTIONALS;i++)
    fun->settings[i] = xcint_param_lookup((enum xc_parameter)i)->default_value;
  fun->nr_active_functionals = 0;
  fun->eval_fun = 0;
  return fun;
}

void xc_free_functional(xc_functional fun)
{
  free(fun);
}

void xc_eval_vec(xc_functional fun, int nr_points,
		 const double *density,
		 int density_pitch,
		 double *result,
		 int result_pitch)
{
  if (fun->mode == XC_MODE_UNSET)
    xcint_die("xc_eval_vec() called before a mode was successfully set",0);
  if (fun->vars == XC_VARS_UNSET)
    xcint_die("xc_eval_vec() called before variables were successfully set",0);
  if (fun->order == -1 && fun->mode != XC_POTENTIAL)
    xcint_die("xc_eval_vec() called before the order was successfully set",0);
  if (fun->eval_fun == 0)
    xcint_die("xc_eval_vec() internal error: no evaluator!",0);
  for (int i=0;i<nr_points;i++)
    fun->eval_fun(fun,density+i*density_pitch,result+i*result_pitch);
}

void xc_eval(xc_functional fun,
	     const double *density,
	     double *result)
{
  if (fun->mode == XC_MODE_UNSET)
    xcint_die("xc_eval() called before a mode was successfully set",0);
  if (fun->vars == XC_VARS_UNSET)
    xcint_die("xc_eval() called before variables were successfully set",0);
  if (fun->order == -1 && fun->mode != XC_POTENTIAL)
    xcint_die("xc_eval() called before the order was successfully set",0);
  if (fun->eval_fun == 0)
    xcint_die("xc_eval() internal error: no evaluator!",0);
  fun->eval_fun(fun,density,result);
}

void xc_set_mode(xc_functional fun, enum xc_mode mode)
{
  fun->mode = mode;
}

void xc_set(xc_functional fun, int item, double value)
{
  if (item >= 0 && item < XC_NR_FUNCTIONALS)
    {
      fun->settings[item] = value;
      fun->active_functionals[fun->nr_active_functionals++] = (enum xc_functional_id)item;
#warning Manage active functional list properly in xc_set  
    }
  else if (item < XC_NR_PARAMETERS_AND_FUNCTIONALS)
    {
      fun->settings[item] = value;
    }
  else
    {
      xcint_die("Invalid item to xc_set():",item);
    }
}

int xc_output_length(xc_functional fun)
{
  if (fun->mode == XC_MODE_UNSET)
    xcint_die("xc_output_length() called before a mode was succesfully set",0);
  if (fun->vars == XC_VARS_UNSET)
    xcint_die("xc_output_length() called before variables were succesfully set",0);
  if (fun->order == -1)
    xcint_die("xc_output_length() called before the order were succesfully set",0);
  if (fun->mode == XC_PARTIAL_DERIVATIVES)
    {
      return taylorlen(xcint_vars_lookup(fun->vars)->partial_vars,fun->order);
    }
  else if (fun->mode == XC_POTENTIAL)
    {
      if (fun->vars == XC_A || fun->vars == XC_A_GAA_TAUA)
	return 2; // Energy+potential
      else
	return 3; // Spin-resolved potential
    }
  else
    {
      xcint_die("XC_CONTRACTED not implemented in xc_output_length()",0);
      return 0;
    }
}

int xcfun_test(void)
{
  int nfail = 0;
  for (int f=0;f<XC_NR_FUNCTIONALS;f++)
    {
      xc_functional fun = xc_new_functional();
      xc_set(fun,f,1.0);
      const functional_data *fd = xcint_functional_lookup((enum xc_functional_id)f);
      if (fd->test_mode != XC_MODE_UNSET)
	{
	  xc_set_mode(fun, fd->test_mode);
	  if (xc_try_vars(fun, fd->test_vars))
	    {
	      if (xc_try_order(fun,fd->test_order))
		{
		  int n = xc_output_length(fun);
		  double *out = 
		    reinterpret_cast<double *>( malloc(sizeof(*out)*n) );
		  if (!fd->test_in)
		    xcint_die("Functional has no test input!",f);
		  xc_eval(fun,fd->test_in,out);
		  int nerr = 0;
		  for (int i=0;i<n;i++)
		    if (fabs(out[i] - fd->test_out[i]) > 
			fabs(fd->test_out[i]*fd->test_threshold))
		      nerr++;
		  if (nerr > 0)
		    {
		      fprintf(stderr,"Error detected in functional %s with tolerance %g:\n",
			      fd->name, fd->test_threshold);
		      fprintf(stderr,"Abs.Error \tComputed              Reference\n");
		      for (int i=0;i<n;i++)
			{
			  fprintf(stderr,"%.1e",fabs(out[i]-fd->test_out[i]));
			  fprintf(stderr,"    %+.16e \t%+.16e",out[i],fd->test_out[i]);
			  if (fabs(out[i] - fd->test_out[i]) > 
			      fabs(fd->test_out[i]*fd->test_threshold))
			    fprintf(stderr," *");
			  fprintf(stderr,"\n");
			}
		      nfail++;
		    }
		  free((void *)out);
		}
	      else
		{
		  fprintf(stderr,"Functional %i not supporting its own test (order)\n",f);
		  nfail++;
		}
	    }
	  else
	    {
	      fprintf(stderr,"Functional %i not supporting its own test (vars)\n",f);
	      nfail++;
	    }
	}
      xc_free_functional(fun);
    }
  return nfail;
}
