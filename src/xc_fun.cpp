#include <cstdlib>
#include <cstring>
#ifndef NO_STDCXX
#include <iostream>
using namespace std;
#endif

#include "xc_fun.h"
#include "taylor.h"
#include "misc.h"
#include "constants.h"

#include "densvars.h"

// Each functional implements a specialization of this template.
template<int functional_id>
struct functional
{
  // Return the name used when selecting this functional.
  static const char *get_name(void);
  // Return a free-form description of the reference and also
  // of any assumptions made in the implementation.
  static const char *get_reference(void);
  // Return the XC correlation energy density,
  // in atomic units (Eh/b^3).
  template<class num>
  static num energy(const densvars<num> &d);
  // Returns 0 if test is ok, otherwise an error code.
  // In particular return -1 if no test is available.
  static int test(void); 
  // Which modes to compile for this functional
  // LDA's should also support GGA modes, but not v.v.
  enum { supported_modes = XC_ALL_GGA }; 
  // Maximum order for this functional. Typically XC_MAX_ORDER
  // but can be lower if you want to save compilation time.
  enum { max_order = XC_MAX_ORDER };
};

#include "testing.h"
#include "settings.h"

#include "lda_fun.h"
#include "blyp_fun.h"
#include "pw92_fun.h"
#include "pbe_fun.h"
#include "lda_sr.h"
#include "tf_kin_fun.h"
#include "pw91_fun.h"

const char *xc_splash(void)
{
  return "XCFun library version "
    XC_FUN_VERSION
    ", Copyright (C) 2009-2010\n"
    " Ulf Ekstrom <uekstrom@gmail.com> and contributors.\n"
    " See http://admol.org/xcfun for a full list of authors.\n\n"
    " This is free software; see the source code for copying conditions.\n"
    " There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or\n"
    " FITNESS FOR A PARTICULAR PURPOSE.\n";
}

// Set up the tables only if we are building the library,
// because the error messages makes developing functionals
// nigh impossible.
#ifdef  XCFUN_LIB
#include "dispatch.h"

static int fun_nvars(enum xc_mode mode)
{
  switch(mode)
    {
    case XC_ALDA: // Fully spin polarized LDA
    case XC_RLDA: // Unpolarized LDA
      return 1;
    case XC_ABLDA: // LDA with alpha/beta densities
    case XC_RSLDA: // LDA with total density/polarization density
    case XC_AGGA:  // Fully spin polarized GGA
    case XC_RGGA:  // Unpolarized GGA
      return 2;
    case XC_ABGGA: // Generic GGA alpha,beta variables 
    case XC_RSGGA: // Generic GGA R,S variables
      return 5;
    case XC_ABMGGA: // Alpha/beta meta gga
    case XC_RSMGGA:
      return 7;
    default:
      assert(0 && "Invalid xc_mode");
      return -1;
    }
}

static int taylor_len(int nvar, int ndeg)
{
  int len = 1;
  for (int k=1;k<=nvar;k++)
    {
      len *= ndeg + k;
      len /= k;
    }
  return len;
}

int xc_len(enum xc_mode mode, int order)
{
  return taylor_len(fun_nvars(mode),order);
}


int xc_index(enum xc_mode mode, const int exponents[])
{
  int nvar = fun_nvars(mode);
  int N = 0;
  for (int i=0;i<nvar;i++)
    N += exponents[i];
  int i = 0, idx = 0;
  N--;
  while (N >= 0)
    {
      idx += taylor_len(nvar-i,N);
      N -= exponents[i];
      i++;
    }
  return idx;
}

//#include "fortran_interface.h"
#endif

template<int fun_id>
static int xc_parse_helper(const char *functional_name)
{
  if (strcmp(functional_name,functional<fun_id>::get_name()) == 0)
    return fun_id;
  else
    return xc_parse_helper<fun_id+1>(functional_name);
}

template<>
int xc_parse_helper<XC_FUNLIST_LEN>(const char *functional_name)
{
  return -1;
}

// Try to get the id number of the functional with
// functional_name. Return -1 if no functional was found.
int xc_fun_id_by_name(const char *functional_name)
{
  return xc_parse_helper<0>(functional_name);
}

static int fun_nvars(enum xc_mode mode)
{
  switch(mode)
    {
    case XC_A:
    case XC_R:
      return 1;
    case XC_AB:
    case XC_RS:
    case XC_AF:
    case XC_RZ:
      return 2;
    case XC_ABFGH:
    case XC_RSZXY:
      return 5;
    default:
      assert(0 && "Invalid xc_mode");
      return -1;
    }
}			    

#if 0


template<class num, int N>
static void fun_helper(num *result, const num d[])
{
  if (fun_nvars(current_mode) == 1)
    {
      taylor<num,1,N> &res = *reinterpret_cast<taylor<num,1,N> *>(result);
      densvars<taylor<num,1,N> > dv;
      if (current_mode == XC_R)
	setdens_r(dv,d);
      sum_functionals(res,dv);
    }
  else if (fun_nvars(current_mode) == 2)
    {
      taylor<num,2,N> &res = *reinterpret_cast<taylor<num,2,N> *>(result);
      densvars<taylor<num,2,N> > dv;
      switch(current_mode)
	{
	case XC_RS: setdens_rs(dv,d); break;
	case XC_AB: setdens_ab(dv,d); break;
	case XC_RZ: setdens_rz(dv,d); break;
	default:
	  break;
	}
      sum_functionals(res,dv);
    }
  else if (fun_nvars(current_mode) == 5)
    {
      taylor<num,5,N> &res = *reinterpret_cast<taylor<num,5,N> *>(result);
      densvars<taylor<num,5,N> > dv;
      if (current_mode == XC_RSZXY)
	setdens_rszxy(dv,d);
      else
	setdens_abfgh(dv,d);
      sum_functionals(res,dv);
#ifdef XC_PRECONDITION
      // Reciprocal scaling factors
      taylor<num,5,N> rscale;
      num scales[5];
      if (current_mode == XC_RSZXY)
	{
	  scales[0] = 1/d[0];
	  scales[1] = 1/d[0];
	  scales[2] = 1/d[2];
	  scales[3] = 1/d[2];
	  scales[4] = 1/d[2];
	}
      else if (current_mode == XC_ABFGH)
	{
	  scales[0] = 1/d[0];
	  scales[1] = 1/d[1];
	  scales[2] = 1/d[2];
	  scales[3] = 1/d[3];
	  scales[4] = 2/(d[2]+d[3]);
	}
      polyterms(rscale,scales);
      for (int i=0;i<res.size;i++)
	res[i] *= rscale[i];
#endif
    }
  else
    {
      assert(0 && "BUG: Invalid functional type");
    }
}

void xc_eval(double *result, int order, const double densvars[])
{
  assert(order>=0);
  if (order > XC_MAX_ORDER)
    {
      cerr << "Error, xc_eval called with order > XC_MAX_ORDER" << endl;
      cerr << "Edit xc_fun.h and recompile." << endl;
      return;
    }
  helper_tab[order](result,densvars);
}


#ifdef WITH_QD
void xc_eval_qd(qd_real *result, int order, const qd_real densvars[])
{
  assert(order>=0);
  assert(order<=XC_MAX_ORDER);
  helper_tab_qd[order](result,densvars);
}
#endif



#endif




#ifndef XCFUN_LIB

template<int functional_id>
void test_functional(void)
{
  assert(xc_fun_id_by_name(functional<functional_id>::get_name()) == 
	 functional_id);
  cout << "   " << functional<functional_id>::get_name() << " ";
  int res = functional<functional_id>::test();
  if (res == 0)
    cout << "\tOk" << endl;
  else if (res == -1)
    cout << "\tNo test" << endl;
  else
    cout << "\tFailed (error " << res << ")" << endl;
  test_functional<functional_id+1>();
}

template<>
void test_functional<XC_FUNLIST_LEN>(void)
{
}

int main(void)
{
  cout << xc_splash() << endl;
  cout << "Testing functionals:" << endl;
  test_functional<0>();
}

#endif
