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
  static const char *get_name(void);
  static const char *get_reference(void);
  template<class num>
  static num energy(const densvars<num> &d);
  // Returns 0 if test is ok, otherwise an error code.
  // In particular return -1 if no test is available.
  static int test(void); 
  // Which modes to compile for this functional
  enum { supported_modes = XC_ALL_GGA }; 
  // Maximum order for this functional
  enum { max_order = XC_MAX_ORDER }; 
};

#include "testing.h"
#include "settings.h"

#include "lda_fun.h"
#include "blyp_fun.h"
#include "pw92_fun.h"
#include "pbe_fun.h"
#include "lda_sr.h"
/*#include "tf_kin_fun.h"
#include "pw91_fun.h"
*/

// _adds_ the functional results weighted with weight.

typedef void (*wrapped_functional)(void *result, const void *dv, const void *weight);

// Dispatch table for all supported functional types and orders
static wrapped_functional double_tab[XC_FUNLIST_LEN][XC_MODES_LEN][XC_MAX_ORDER+1] = {{{0}}};
static const char *fun_name_tab[XC_FUNLIST_LEN] = {0};//"knas","boll"}; 

#ifdef  XCFUN_LIB
#include "dispatch.h"
#endif

#if 0



template<class scalar, enum xc_mode mode, int Ndeg>
class composite_functional
{
public:
  composite_functional(void) : term(0), next(0), weight(0) {}
  int add_part(enum xc_functionals fun, scalar weight)
  {
    if (m_term == 0)
      {
	m_term = functional_picker<scalar,mode,Ndeg,0>::get_by_id(fun);
	if (not m_term)
	  return -1;
	m_weight = weight;
      }
    else 
      {
	if (not m_next)
	  m_next = new composite_functional<scalar,mode,Ndeg>;
	m_next.add_part(fun,weight);
      }
    return 0;
  }
  int add_part(const char *name, scalar weight)
  {
    enum xc_functionals id = functional_picker<scalar,mode,Ndeg,0>::get_id_by_name(name);
    if (id < 0)
      return -1;
    add_part(id,weight);
  }
private:
  typedef taylor<scalar,xc_mode_info<mode>::Nvar,Ndeg> ttype;
  typedef ttype (*ftype)(const densvars<ttype> &);

  scalar m_weight;
  ttype (*m_term)(const densvars<ttype> &d);
  composite_functional<scalar,mode,Ndeg> *m_next;
};

#endif
#if 0
//Did we already set up tables?
static int is_setup = 0;
//Which primary variables will we use?
static enum xc_mode current_mode = XC_A;

static void setup(void);

// Use strcmp to not depend on c++ standard library
int xc_parse_functional(const char *functional_name)
{
  for (int i=0;i<XC_PARAMS_LEN;i++)
    if (strcmp(functional_name,functional_names[i]) == 0)
      return i;
  return -1;
}

void xc_set_functional(enum xc_mode mode, 
		       const double params[XC_PARAMS_LEN])
{
  if (!is_setup)
    setup();
  current_mode = mode;
  for (int i=0;i<XC_PARAMS_LEN;i++)
    xc_params[i] = params[i];
}

extern "C"
void xc_set_functional_fortran_(const int *mode, const double *params, int *len)
{
  assert(*len == XC_PARAMS_LEN && "C++/Fortran headers out of sync");
  xc_set_functional((enum xc_mode)(*mode),params);
}

#ifndef NO_STDCXX
extern "C"
void fortran_report_xc_functional_()
{
  // this is useful for debugging mpi
  // it will report the set functional mix
  // which may not be set on a mpi worker
  for (int i = 0; i < XC_PARAMS_LEN; i++)
    cout << "functional portion:" << i << " weight:" << xc_params[i] << endl;
}
#endif



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

#ifdef WITH_SINGLE
void xc_eval_single(float *result, int order, const float densvars[])
{
  assert(order>=0);
  assert(order<=XC_MAX_ORDER);
  helper_tab_single[order](result,densvars);
}
#endif

#ifdef WITH_QD
void xc_eval_qd(qd_real *result, int order, const qd_real densvars[])
{
  assert(order>=0);
  assert(order<=XC_MAX_ORDER);
  helper_tab_qd[order](result,densvars);
}
#endif

extern "C"
void xc_eval_fortran_(double *result, int *order, const double *densvars)
{
  xc_eval(result,*order,densvars);
}

#ifdef WITH_QD
extern "C"
void xc_eval_fortran_qd_(double *result, int *order, const double *densvars)
{
  qd_real qd_result[xc_len(current_mode,*order)];
  qd_real qd_densvars[fun_nvars(current_mode)];
  for (int i=0;i<fun_nvars(current_mode);i++)
    qd_densvars[i] = densvars[i];
  xc_eval_qd(qd_result,*order,qd_densvars);
  for (int i=0;i<xc_len(current_mode,*order);i++)
    result[i] = to_double(qd_result[i]);
}
#endif


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

extern "C"
int xc_len_fortran_(int *mode, int *order)
{
  return xc_len((enum xc_mode)*mode,*order);
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

// Fortran version
extern "C"
int xc_index_fortran_(const int *mode, const int exponents[])
{
  return 1+xc_index((enum xc_mode)*mode,exponents);
}

template<int N>
void setup_order(void)
{
  helper_tab[N] = fun_helper<double,N>;
#ifdef WITH_QD
  helper_tab_qd[N] = fun_helper<qd_real,N>;
#endif
#ifdef WITH_SINGLE
  helper_tab_single[N] = fun_helper<float,N>;
#endif
  setup_order<N-1>();
}

#endif

#ifndef XCFUN_LIB

template<int functional_id>
void test_functional(void)
{
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
  //setup_tables();
  /*
  cout << "Available functionals (orders):" << endl;
  for (int i=0;i<XC_FUNLIST_LEN;i++)
    if (fun_name_tab[i])
      {
	cout << "   " << fun_name_tab[i] << " \t";
	for (int j=0;j<=XC_MAX_ORDER;j++)
	  if (double_tab[i][bitpos<XC_ABGGA>::pos][j])
	    cout << " " << j;
	cout << endl;
      }
  */
  cout << "Testing functionals:" << endl;
  test_functional<0>();
}

#endif
