#include "xcint.h"
#include <cstdlib>

template<class T, int F>
struct spec_table
{
  static void setup(T (**table)(const densvars<T> &d))
  {
    table[F] = erg<F,T>::f;
    //    printf("table[%i] = %p\n",F,table[F]);
    spec_table<T,F-1>::setup(table);
  }
};

template<class T>
struct spec_table<T,-1>
{
  static void setup(T (**table)(const densvars<T> &d))
  {
  }
};

// Return a pointer to a specialization for the functional f on the type T,
// or NULL if no such specialization exists.
template<class T>
T (*xcint_spec(enum xc_functional_id f))(const densvars<T> &d)
{
  static T (**table)(const densvars<T> &d) = 0;
  if (!table)
    {
      // TODO: NOT THREADSAFE
      table = 
	reinterpret_cast<T (**)(const densvars<T> &d)>
	(malloc(sizeof(T (*)(const densvars<T> &d))*XC_NR_FUNCTIONALS));
      spec_table<T,XC_NR_FUNCTIONALS-1>::setup(table);
    }
  return table[f];
}

template<int MODE, int VARS, int ORDER>
struct evaluator_instance
{
  static void eval_fun(xc_functional_obj *f, const double *input, double *output);
  static evaluator get_evaluator(enum xc_mode mode,
				 enum xc_vars vars,
				 int order)
  {
    if (mode == MODE && vars == VARS && order == ORDER)
      return evaluator_instance<MODE,VARS,ORDER>::eval_fun;
    else
      return evaluator_instance<MODE,VARS,ORDER-1>::get_evaluator(mode,vars,order);
  }
};



evaluator xcint_get_evaluator(enum xc_mode mode,
			      enum xc_vars vars,
			      int order)
{
  return evaluator_instance<XC_NR_MODES-1,XC_NR_VARS-1,XC_MAX_ORDER>::get_evaluator(mode,vars,order);
}

/* --- --- --- loops below --- --- --- */


template<int MODE, int VARS>
struct evaluator_instance<MODE,VARS,-1>
{
  static evaluator get_evaluator(enum xc_mode mode, enum xc_vars vars, int order)
  {
    return evaluator_instance<MODE,VARS-1,XC_MAX_ORDER>::get_evaluator(mode,vars,order);
  }
};

template<int MODE, int ORDER>
struct evaluator_instance<MODE,-1,ORDER>
{
  static evaluator get_evaluator(enum xc_mode mode, enum xc_vars vars, int order)
  {
    return evaluator_instance<MODE-1,XC_NR_VARS-1,XC_MAX_ORDER>::get_evaluator(mode,vars,order);
  }
};

template<int MODE>
struct evaluator_instance<MODE,-1,-1>
{
  static evaluator get_evaluator(enum xc_mode mode, enum xc_vars vars, int order)
  {
    return evaluator_instance<MODE-1,XC_NR_VARS-1,XC_MAX_ORDER>::get_evaluator(mode,vars,order);
  }
};

template<int VARS, int ORDER>
struct evaluator_instance<-1,VARS,ORDER>
{
  static evaluator get_evaluator(enum xc_mode mode, enum xc_vars vars, int order)
  {
    return 0;
  }
};

template<>
struct evaluator_instance<-1,-1,-1>
{
  static evaluator get_evaluator(enum xc_mode mode, enum xc_vars vars, int order)
  {
    return 0;
  }
};

/* --- --- --- ---       --- --- --- */

// This is a replacement for copysign from C99
static inline ireal_t cpsign(ireal_t x, ireal_t y)
{
  if (x<0)
    {
      if (y<0)
	return x;
      else
	return -x;
    }
  else
    {
      if (y<0)
	return -x;
      else
	return x;
    }
}

// Here VARS are the variables that are actually used during the functional evaluation, i.e.
// not laplacian for gga densities
template<int VARS, class T>
void xcint_setup_vars(densvars<T> &dv, const double *d)
{
  if (VARS == XC_A || VARS == XC_A_GAA || VARS == XC_A_GAA_TAUA)
    {
#ifdef XC_NO_REGULARIZATION
      dv.a = T(d[0],0);
#else
      dv.a = T(d[0] > XC_TINY_DENS ? d[0] : XC_TINY_DENS,0);
#endif
      dv.b = 0;
      dv.n = dv.a;
      dv.s = dv.a;
      dv.zeta = 1;
      dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
      dv.n_m13 = pow(dv.n,-1.0/3.0);
      dv.a_43 = pow(dv.a,4.0/3.0);
      dv.b_43 = 0;
      if (VARS != XC_A)
	{
#ifdef XC_NO_REGULARIZATION
	  dv.gaa = T(d[1],1);
#else
	  dv.gaa = T(d[1] >= 0 ? d[2] : 0,2);
#endif
	  dv.gab = 0;
	  dv.gbb = 0;  
	  dv.gnn  = dv.gaa;
	  dv.gss  = dv.gaa;
	  dv.gns  = dv.gaa;
	  if (VARS == XC_A_GAA_TAUA)
	    {
#ifdef XC_NO_REGULARIZATION
	      dv.taua = T(d[2], 2);
#else
	      dv.taua = T(d[2] >= XC_TINY_DENS ? d[2] : XC_TINY_DENS, 2);
#endif
	      dv.taub = 0;
	      dv.tau = dv.taua + dv.taub;
	    }
	}
    }
  if (VARS == XC_A_B || VARS == XC_A_B_GAA_GAB_GBB || VARS == XC_A_B_GAA_GAB_GBB_TAUA_TAUB)
    {
#ifdef XC_NO_REGULARIZATION
      dv.a = T(d[0],0);
      dv.b = T(d[1],1);
#else
      dv.a = T(d[0] > XC_TINY_DENS ? d[0] : XC_TINY_DENS,0);
      dv.b = T(d[1] > XC_TINY_DENS ? d[1] : XC_TINY_DENS,1);
#endif  
      dv.n = dv.a+dv.b;
      dv.s = dv.a-dv.b;
      dv.zeta = dv.s/dv.n;
      dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
      dv.n_m13 = pow(dv.n,-1.0/3.0);
      dv.a_43 = pow(dv.a,4.0/3.0);
      dv.b_43 = pow(dv.b,4.0/3.0);    
      if (VARS != XC_A_B)
	{
#ifdef XC_NO_REGULARIZATION
	  dv.gaa = T(d[2],2);
	  dv.gab = T(d[3],3);
	  dv.gbb = T(d[4],4);
#else
	  dv.gaa = T(d[2] >= 0 ? d[2] : 0,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
	  dv.gab = T(d[3],3); 
#else
	  dv.gab = T(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
	  dv.gbb = T(d[4] >= 0 ? d[4] : 0,4);
#endif  
	  dv.gnn  = dv.gaa + 2*dv.gab + dv.gbb; 
	  dv.gss  = dv.gaa - 2*dv.gab + dv.gbb;
	  dv.gns  = dv.gaa - dv.gbb;
	  if (VARS == XC_A_B_GAA_GAB_GBB_TAUA_TAUB)
	    {
#ifdef XC_NO_REGULARIZATION
	      dv.taua = T(d[5], 5);
	      dv.taub = T(d[6], 6);
#else
	      dv.taua = T(d[5] >= XC_TINY_DENS ? d[5] : XC_TINY_DENS, 5);
	      dv.taub = T(d[6] >= XC_TINY_DENS ? d[6] : XC_TINY_DENS, 6);
#endif
	      dv.tau = dv.taua + dv.taub;
	    }
	}
    }
  else
    {
      xcint_die("VARS value not implemented",VARS);
    }
}

template<int MODE, int VARS, int ORDER>
void evaluator_instance<MODE,VARS,ORDER>::eval_fun(xc_functional_obj *f, const double *input, double *output)
{
  if (MODE == XC_PARTIAL_DERIVATIVES)
    {
      typedef taylor<ireal_t,vars_info<VARS>::nr_variables,ORDER> ttype;
      ttype *result = reinterpret_cast<ttype *>(output);
      densvars<ttype> d;
      xcint_setup_vars<VARS,ttype>(d,input);
      d.parent = f;
      *result = 0;
      for (int i=0;i<f->nr_active_functionals;i++)
	{
	  if (xcint_spec<ttype>(f->active_functionals[i]))
	    *result += f->settings[f->active_functionals[i]] *
	      xcint_spec<ttype>(f->active_functionals[i])(d);
	  else
	    xcint_die("Functional does not have requested specialization",f->active_functionals[i]);
	}
      result->deriv_facs();
    }
  else if (MODE == XC_POTENTIAL)
    {
      // First derivatives needed for lda, second for gga
      // TODO: select proper type here
      typedef taylor<ireal_t,vars_info<VARS>::nr_variables,vars_info<VARS>::pot_order> ttype;
      ttype out = 0;
      densvars<ttype> d;
      xcint_setup_vars<VARS,ttype>(d,input);
      d.parent = f;
      for (int i=0;i<f->nr_active_functionals;i++)
	{
	  if (xcint_spec<ttype>(f->active_functionals[i]))
	    out += f->settings[f->active_functionals[i]] *
	      xcint_spec<ttype>(f->active_functionals[i])(d);
	  else
	    xcint_die("Functional does not have requested specialization",f->active_functionals[i]);
	}
      output[0] = out[0];
      if (VARS == XC_A or VARS == XC_N)
	{
	  output[1] = out[1];
	}
      else if (VARS == XC_A_B or VARS == XC_N_S)
	{
	  output[1] = out[1];
	  output[2] = out[2];
	}
      else if (VARS == XC_A_B_GAA_GAB_GBB)
	{
	  const int gaa = 2, gab = 3, gbb = 4, lapa = 5, lapb = 6;
	  output[1] = out[1];
	  output[2] = out[2];

	  output[1]  = out[XC_D10000];
	  output[1] -= 2*input[lapa]*out[XC_D00100] + input[lapb]*out[XC_D00010];
	  output[1] -= 2*(out[XC_D10100]*input[gaa]   + 
			out[XC_D01100]*input[gab] +
			out[XC_D00200]*(2*input[lapa]*input[gaa]) +
			out[XC_D00110]*(input[lapa]*input[gab] + input[lapb]*input[gaa]) +
			out[XC_D00101]*(2*input[lapb]*input[gab]) 
			); 
	  output[1] -= (out[XC_D10010]*input[gab] +
		      out[XC_D01010]*input[gbb] +
		      out[XC_D00110]*(2*input[lapa]*input[gab]) +
		      out[XC_D00020]*(input[lapb]*input[gab] + input[lapa]*input[gbb]) +
		      out[XC_D00011]*(2*input[lapb]*input[gbb])); 

	  output[2]  = out[XC_D01000];
	  output[2] -= 2*input[lapb]*out[XC_D00001] + input[lapa]*out[XC_D00010];
	  output[2] -= 2*(out[XC_D01001]*input[gbb]   + 
			out[XC_D10001]*input[gab]  +
			out[XC_D00002]*(2*input[lapb]*input[gbb]) +
			out[XC_D00011]*(input[lapb]*input[gab] + input[lapa]*input[gbb]) +
			out[XC_D00101]*(2*input[lapa]*input[gab])  ); 
	  output[2] -= (out[XC_D01010]*input[gab] +
			out[XC_D10010]*input[gaa] +
			out[XC_D00011]*(2*input[lapb]*input[gab]) +
			out[XC_D00020]*(input[lapa]*input[gab] + input[lapb]*input[gaa]) +
			out[XC_D00110]*(2*input[lapa]*input[gaa])); 
	}
      else
	{
	  xcint_die("XC_POTENTIAL not implemented for this functional and variables",0);
	}
    }
  else if (MODE == XC_CONTRACTED)
    {
      xcint_die("XC_CONTRACT not implemented",0);
    }
  else
    {
      xcint_die("Invalid MODE in evaluator_instance",MODE);
    }
}

// Assuming fun is set up with a valid set of (mode, vars, order)
// assign an evaluator to use
static void xcint_pick_evaluator(xc_functional_obj *fun)
{
  fun->eval_fun = evaluator_instance<XC_NR_MODES-1,XC_NR_VARS-1,XC_MAX_ORDER>
    ::get_evaluator(fun->mode,fun->vars,fun->order);
  //  printf("Evaluator set to %p\n",fun->eval_fun);
}

template<int VARS, int ORDER>
bool xcint_try_helper(xc_functional fun, enum xc_vars vars, int order)
{
  if (vars == VARS && order == ORDER)
    {
      for (int i=0;i<fun->nr_active_functionals;i++)
	if (xcint_spec<taylor<ireal_t,vars_info<VARS>::nr_variables,ORDER> >(fun->active_functionals[i]) == 0)
	  {
	    printf("looked for implementation <ireal_t, %i, %i> for functional %i, found nothing.\n",
		   vars_info<VARS>::nr_variables,ORDER, fun->active_functionals[i]);
	    return false;
	  }
      return true;
    }
  else
    {
      return xcint_try_helper
	<ORDER == 0 ? VARS-1 : VARS, ORDER == 0 ? XC_MAX_ORDER : ORDER - 1>(fun,vars,order);
    }
}

template<>
bool xcint_try_helper<-1,XC_MAX_ORDER>(xc_functional fun, enum xc_vars vars, int order) { return false; }


int xc_try_vars(xc_functional fun, enum xc_vars vars)
{
  int order = fun->order;
  if (order < 0)
    order = 0;
  if (!xcint_try_helper<XC_NR_VARS-1,XC_MAX_ORDER>(fun, vars, order))
    return 0;
  fun->vars = vars;
  if (fun->mode != XC_MODE_UNSET && fun->order >= 0)
    xcint_pick_evaluator(fun);
  return 1;
}


int xc_try_order(xc_functional fun, int order)
{
  enum xc_vars vars = fun->vars;
  if (vars == XC_VARS_UNSET)
    vars = XC_A;
  if (!xcint_try_helper<XC_NR_VARS-1,XC_MAX_ORDER>(fun, vars, order))
    return 0;
  fun->order = order;
  if (fun->mode != XC_MODE_UNSET && fun->vars != XC_VARS_UNSET && fun->order >= 0)
    xcint_pick_evaluator(fun);
  return 1;
}
