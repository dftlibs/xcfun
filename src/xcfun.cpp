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
