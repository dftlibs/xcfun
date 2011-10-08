#include <cstdio>
#include <cstdlib>
#include "xcint.hpp"

xc_functional xc_new_functional_not_macro(int checksum)
{
  xcint_assure_setup();
  if (checksum != XCFUN_CHECKSUM)
    xcint_die("Header/library inconsistency detected, aborting",checksum);
  xc_functional fun = (xc_functional)malloc(sizeof*fun);
  if (!fun)
    xcint_die("Out of memory in xc_new_functional()",0);
  fun->mode = XC_MODE_UNSET;
  fun->vars = XC_VARS_UNSET;
  fun->order = -1;
  fun->depends = 0;
  for (int i=0;i<XC_NR_FUNCTIONALS;i++)
    fun->settings[i] = 0;
  for (int i=XC_NR_FUNCTIONALS;i<XC_NR_PARAMETERS_AND_FUNCTIONALS;i++)
    fun->settings[i] = xcint_params[i].default_value;
  fun->nr_active_functionals = 0;
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
  for (int i=0;i<nr_points;i++)
    xc_eval(fun,density+i*density_pitch,result+i*result_pitch);
}

void xc_set(xc_functional fun, int item, double value)
{
  if (item >= 0 && item < XC_NR_FUNCTIONALS)
    {
      fun->settings[item] = value;
      fun->active_functionals[fun->nr_active_functionals++] = &xcint_funs[item];
      fun->depends |= xcint_funs[item].depends;
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

double xc_get(xc_functional fun, int item)
{
  if (item >= 0 && item < XC_NR_PARAMETERS_AND_FUNCTIONALS)
    {
      return fun->settings[item];
    }
  else
    {
      xcint_die("Invalid item to xc_get():",item);
      return 0;
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
      return taylorlen(xcint_vars[fun->vars].len,fun->order);
    }
  else if (fun->mode == XC_POTENTIAL)
    {
      if (fun->vars == XC_A || fun->vars == XC_A_2ND_TAYLOR)
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
  int nfail = 0, res;
  for (int f=0;f<XC_NR_FUNCTIONALS;f++)
    {
      xc_functional fun = xc_new_functional();
      xc_set(fun,f,1.0);
      const functional_data *fd = &xcint_funs[f];
      if (fd->test_mode != XC_MODE_UNSET)
	{
	  if ((res = xc_eval_setup(fun,fd->test_vars,fd->test_mode,fd->test_order)) == 0)
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
			  fd->symbol, fd->test_threshold);
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
	      fprintf(stderr,"Functional %s not supporting its own test, error %i\n",fd->symbol,res);
	      nfail++;
	    }
	}
      else
	{
	  fprintf(stderr,"%s has no test\n",fd->symbol);
	}
      xc_free_functional(fun);
    }
  return nfail;
}

double xcfun_version(void)
{
  return 1.99;
}

const char *xcfun_splash(void)
{
  return 
    "XCFun DFT library Copyright 2009-2011 Ulf Ekstrom and contributors.\n"
    "See http://admol.org/xcfun for more information.\n\n"
    "This is free software; see the source code for copying conditions.\n"
    "There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or\n"
    "FITNESS FOR A PARTICULAR PURPOSE. For details see the documentation.\n"
    "Scientific users of this library should cite\n"
    "U. Ekstrom, L. Visscher, R. Bast, A. J. Thorvaldsen and K. Ruud;\n"
    "J.Chem.Theor.Comp. 2010, DOI: 10.1021/ct100117s\n";
}

void xc_eval(xc_functional_obj *f, const double *input, double *output)
{
  if (f->mode == XC_MODE_UNSET)
    xcint_die("xc_eval() called before a mode was successfully set",0);
  if (f->vars == XC_VARS_UNSET)
    xcint_die("xc_eval() called before variables were successfully set",0);
  if (f->order == -1 && f->mode != XC_POTENTIAL)
    xcint_die("xc_eval() called before the order was successfully set",0);
  if (f->mode == XC_PARTIAL_DERIVATIVES)
    {
      switch (f->order)
	{
	case 0:
	  {
	    typedef ctaylor<ireal_t,0> ttype;
	    int inlen = xcint_vars[f->vars].len;
	    ttype in[inlen], out = 0;
	    for (int i=0;i<inlen;i++)
	      in[i] = input[i];
	    densvars<ttype> d(f,in);
	    for (int i=0;i<f->nr_active_functionals;i++)
	      out +=  f->settings[f->active_functionals[i]->id]
		      * f->active_functionals[i]->fp0(d);
	    output[0] = out.get(CNST);
	  }
	  break;
#if XC_MAX_ORDER >= 1
	case 1:
	  {
	  int inlen = xcint_vars[f->vars].len;
	  {
	    typedef ctaylor<ireal_t,2> ttype2;
	    ttype2 in2[inlen], out2 = 0;
	    for (int i=0;i<inlen;i++)
	      in2[i] = input[i];
	    for (int j=0;j<inlen/2;j++)
	      {
		in2[2*j].set(VAR0,1);
		in2[2*j+1].set(VAR1,1);
		densvars<ttype2> d(f,in2);
		out2 = 0;
		for (int i=0;i<f->nr_active_functionals;i++)
		  out2 += f->settings[f->active_functionals[i]->id]
		    * f->active_functionals[i]->fp2(d);
		in2[2*j] = input[2*j];
		in2[2*j+1] = input[2*j+1];
		output[2*j+1] = out2.get(VAR0); // First derivatives
		output[2*j+2] = out2.get(VAR1); // First derivatives
	      }
	    output[0] = out2.get(CNST); // Energy
	  }
	  if (inlen & 1)
	  {
	    typedef ctaylor<ireal_t,1> ttype;
	    int inlen = xcint_vars[f->vars].len;
	    ttype in[inlen], out = 0;
	    for (int i=0;i<inlen;i++)
	      in[i] = input[i];
	    int j = inlen-1;
	      {
		in[j].set(VAR0,1);
		densvars<ttype> d(f,in);
		out = 0;
		for (int i=0;i<f->nr_active_functionals;i++)
		  out += f->settings[f->active_functionals[i]->id]
		    * f->active_functionals[i]->fp1(d);
		in[j] = input[j];
		output[j+1] = out.get(VAR0); // First derivatives
	      }
	      output[0] = out.get(CNST); // Energy
	  }
	  }
	  break;
#if 0
	    typedef ctaylor<ireal_t,1> ttype;
	    int inlen = xcint_vars[f->vars].len;
	    ttype in[inlen], out = 0;
	    for (int i=0;i<inlen;i++)
	      in[i] = input[i];
	    for (int j=0;j<inlen;j++)
	      {
		in[j].set(VAR0,1);
		densvars<ttype> d(f,in);
		out = 0;
		for (int i=0;i<f->nr_active_functionals;i++)
		  out += f->settings[f->active_functionals[i]->id]
		    * f->active_functionals[i]->fp1(d);
		in[j] = input[i];
		output[j+1] = out[1]; // First derivatives
	      }
	    output[0] = out[0]; // Energy
	  }
	  break;
#endif
#endif
#if XC_MAX_ORDER >= 2
	case 2:
	  {
	    typedef ctaylor<ireal_t,2> ttype;
	    int inlen = xcint_vars[f->vars].len;
	    ttype in[inlen], out = 0;
	    for (int i=0;i<inlen;i++)
	      in[i] = input[i];
	    int k = inlen+1;
	    for (int i=0;i<inlen;i++)
	      {
		in[i].set(VAR0,1);
		for (int j=i;j<inlen;j++)
		  {
		    in[j].set(VAR1,1);
		    densvars<ttype> d(f,in);
		    out = 0;
		    for (int n=0;n<f->nr_active_functionals;n++)
		      out += f->settings[f->active_functionals[n]->id]
		      * f->active_functionals[n]->fp2(d); 
		    output[k++] = out.get(VAR0|VAR1); //Second derivative
		    in[j].set(VAR1,0); // slightly pessimized
		  }
		output[i+1] = out.get(VAR0); //First derivative
		in[i] = input[i];
	      }
	    output[0] = out.get(CNST); //Energy
	  }
	  break;
#endif
	default:
	  xcint_die("FIXME: Order too high for partial derivatives in xc_eval",f->order);
	}
    }
  else if (f->mode == XC_CONTRACTED)
    {
#define DOEVAL(N,E) if (f->order == N) {\
      typedef ctaylor<ireal_t,N> ttype; \
      int inlen = xcint_vars[f->vars].len; \
      ttype in[inlen], out = 0; \
      int k = 0; \
      for (int i=0;i<inlen;i++) \
	for (int j=0;j<(1 << f->order);j++) \
	  in[i].set(j,input[k++]);	    \
      densvars<ttype> d(f,in); \
      for (int i=0;i<f->nr_active_functionals;i++) \
	out += f->settings[f->active_functionals[i]->id] \
	  * f->active_functionals[i]->fp##N(d);\
      for (int i=0;i<(1<<f->order);i++) output[i] = out.get(i); } else
      FOR_EACH(XC_MAX_ORDER,DOEVAL,)
      xcint_die("bug! Order too high in XC_CONTRACTED",f->order);
    }
  else if (f->mode == XC_POTENTIAL)
    {
      int inlen = xcint_vars[f->vars].len;
      int npot; // One or two potentials
      if (inlen == 1 || inlen == 10)
	npot = 1;
      else
	npot = 2;

      {
	typedef ctaylor<ireal_t,1> ttype;      
	ttype in[inlen], out = 0;
	for (int i=0;i<inlen;i++)
	  in[i] = input[i];
	for (int j=0;j<npot;j++)
	  {
	    in[j].set(VAR0,1);
	    densvars<ttype> d(f,in);
	    out = 0;
	    for (int i=0;i<f->nr_active_functionals;i++)
	      out += f->settings[f->active_functionals[i]->id]
		* f->active_functionals[i]->fp1(d);
	    in[j] = input[j];
	    output[j+1] = out.get(VAR0); // First derivatives
	  }
	output[0] = out.get(CNST); // Energy
      }

      if (f->depends & XC_GRADIENT) // GGA potential
	{
	  /*
	     v = dE/dn - nabla.dE/dg
	   */
	  typedef ctaylor<ireal_t,2> ttype;
	  ttype in[inlen];
	  if (f->vars == XC_A_2ND_TAYLOR || f->vars == XC_N_2ND_TAYLOR)
	    {
	      // n gx gy gz xx xy xz yy yz zz
	      // 0 1  2  3  4  5  6  7  8  9
	      ttype out = 0;
	      // d/dx
	      in[0] = ttype(input[0],0,input[1]);
	      for (int i=0;i<3;i++)
		in[1+i] = ttype(input[1+i],VAR0,input[4+i]);
	      in[1].set(VAR1,1); // d/dgx
	      {
		densvars<ttype> d(f,in);
		for (int i=0;i<f->nr_active_functionals;i++)
		  out += f->settings[f->active_functionals[i]->id]
		    * f->active_functionals[i]->fp2(d);
	      }

	      // d/dy
	      in[0] = ttype(input[0],VAR0,input[2]);
	      in[1] = ttype(input[1],VAR0,input[5]);
	      in[2] = ttype(input[2],VAR0,input[7]);
	      in[3] = ttype(input[3],VAR0,input[8]);
	      in[2].set(VAR1,1); // d/dgy
	      {
		densvars<ttype> d(f,in);
		for (int i=0;i<f->nr_active_functionals;i++)
		  out += f->settings[f->active_functionals[i]->id]
		    * f->active_functionals[i]->fp2(d);
	      }
	      // d/dz
	      in[0] = ttype(input[0],VAR0,input[3]);
	      in[1] = ttype(input[1],VAR0,input[6]);
	      in[2] = ttype(input[2],VAR0,input[8]);
	      in[3] = ttype(input[3],VAR0,input[9]);
	      in[3].set(VAR1,1); // d/dgz
	      {
		densvars<ttype> d(f,in);
		for (int i=0;i<f->nr_active_functionals;i++)
		  out += f->settings[f->active_functionals[i]->id]
		    * f->active_functionals[i]->fp2(d);
	      }
	      output[1] -= out.get(VAR0|VAR1); // Subtract divergence of dE/dg from lda part of potential
	    }
	  else if (f->vars == XC_A_B_2ND_TAYLOR || f->vars == XC_N_S_2ND_TAYLOR)
	    {
	      xcint_die("AB potential Under construction in xc_eval()",f->mode);
	      ttype out = 0;
	      for (int spin=0;spin<=1;spin++)
		{
		  int spinoff = 10*spin;
		  // d/dx
		  in[0] = ttype(input[0+spinoff],0,input[1+spinoff]);
		  for (int i=0;i<3;i++)
		    in[1+i] = ttype(input[1+i+spinoff],VAR0,input[4+i+spinoff]);
		  in[1].set(VAR1,1); // d/dgx
		  {
		    densvars<ttype> d(f,in);
		    for (int i=0;i<f->nr_active_functionals;i++)
		      out += f->settings[f->active_functionals[i]->id]
			* f->active_functionals[i]->fp2(d);
		  }

		  // d/dy
		  in[0] = ttype(input[0+spinoff],VAR0,input[2+spinoff]);
		  in[1] = ttype(input[1+spinoff],VAR0,input[5+spinoff]);
		  in[2] = ttype(input[2+spinoff],VAR0,input[7+spinoff]);
		  in[3] = ttype(input[3+spinoff],VAR0,input[8+spinoff]);
		  in[2].set(VAR1,1); // d/dgy
		  {
		    densvars<ttype> d(f,in);
		    for (int i=0;i<f->nr_active_functionals;i++)
		      out += f->settings[f->active_functionals[i]->id]
			* f->active_functionals[i]->fp2(d);
		  }
		  // d/dz
		  in[0] = ttype(input[0+spinoff],VAR0,input[3+spinoff]);
		  in[1] = ttype(input[1+spinoff],VAR0,input[6+spinoff]);
		  in[2] = ttype(input[2+spinoff],VAR0,input[8+spinoff]);
		  in[3] = ttype(input[3+spinoff],VAR0,input[9+spinoff]);
		  in[3].set(VAR1,1); // d/dgz
		  {
		    densvars<ttype> d(f,in);
		    for (int i=0;i<f->nr_active_functionals;i++)
		      out += f->settings[f->active_functionals[i]->id]
			* f->active_functionals[i]->fp2(d);
		  }
		  output[1+spin] -= out.get(VAR0|VAR1); // Subtract divergence of dE/dg from lda part of potential
		}
	    }
	}
    }
  else 
    {
      xcint_die("Illegal mode in xc_eval()",f->mode);
    }
}

int xc_eval_setup(xc_functional fun,
		  enum xc_vars vars,
		  enum xc_mode mode,
		  int order)
{
  // Check that vars are enough for the functional
  if ((fun->depends & xcint_vars[vars].provides) != fun->depends)
    {
      printf("depends %i, provides %i\n",fun->depends,xcint_vars[vars].provides);
      return XC_EVARS;
    }
  if ((order < 0 || order > XC_MAX_ORDER) ||
      (mode == XC_PARTIAL_DERIVATIVES && order > 2))
    return XC_EORDER;
  if (mode == XC_POTENTIAL)
    {
      // GGA potential needs full laplacian
      if ((fun->depends & XC_GRADIENT) && !(vars == XC_A_2ND_TAYLOR ||
					    vars == XC_A_B_2ND_TAYLOR ||
					    vars == XC_N_2ND_TAYLOR ||
					    vars == XC_N_S_2ND_TAYLOR))
	{
	  printf("xxx\n");
	  return XC_EVARS | XC_EMODE;
	  
	}
      // No potential for meta gga's
      if (fun->depends & (XC_LAPLACIAN | XC_KINETIC))
	return XC_EMODE;
    }
  fun->mode = mode;
  fun->vars = vars;
  fun->order = order;
  return 0;
}

const char *xc_name(int param)
{
  xcint_assure_setup();
  if (param >= 0 && param < XC_NR_FUNCTIONALS)
    return xcint_funs[param].symbol;
  else if (param < XC_NR_PARAMETERS_AND_FUNCTIONALS)
    return xcint_params[param].symbol;
  else
    return 0;
}

// Describe in one line what the setting does
const char *xc_short_description(int param)
{
  xcint_assure_setup();
  if (param >= 0 && param < XC_NR_FUNCTIONALS)
    return xcint_funs[param].short_description;
  else if (param < XC_NR_PARAMETERS_AND_FUNCTIONALS)
    return xcint_params[param].description;
  else
    return 0;
}
// Long description of the setting, ends with a \n
const char *xc_long_description(int param)
{
  xcint_assure_setup();
  if (param >= 0 && param < XC_NR_FUNCTIONALS)
    return xcint_funs[param].long_description;
  else if (param < XC_NR_PARAMETERS_AND_FUNCTIONALS)
    return xcint_params[param].description;
  else
    return 0;
}
