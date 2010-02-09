#include <cstdlib>
#include <cstring>
#ifndef NO_STDCXX
#include <iostream>
using namespace std;
#endif

#include "taylor.h"

#include "constants.h"

#include "xc_fun.h"
#include "lda_fun.h"
#include "blyp_fun.h"
#include "pbe_fun.h"
#include "tf_kin_fun.h"
#include "pw91_fun.h"

template<class T> 
struct densvars
{
  // Variables for expressing functionals, these are redundant because
  // different functionals have different needs.
  T R, S, Z, X, Y;
  T rhoa, rhob, grada2, gradab, gradb2;
  T zeta; //S/R
};

//Did we already set up tables?
static int is_setup = 0;
//Which primary variables will we use?
static enum xc_mode current_mode = XC_A;
// Functional parameters
static double xc_params[XC_PARAMS_LEN];
//Dispatch table for different orders, set by setup()
static void (*helper_tab[XC_MAX_ORDER+1])(double *, const double *);
#ifdef WITH_QD
static void (*helper_tab_qd[XC_MAX_ORDER+1])(qd_real *, const qd_real *);
#endif
#ifdef WITH_SINGLE
static void (*helper_tab_single[XC_MAX_ORDER+1])(float *, const float *);
#endif
static void setup(void);

/* This must match the order of enum xc_funparams! 
   This is defined as an array of char* to avoid relying
   on the C++ runtime (although the C++ string class
   would be better to use).
*/
static const char *functional_names[] =
  {
    "vwn5c",
    "slaterx",
    "lypc",
    "pw91x",
    "pbec",
    "pbec_dalton",
    "beckex",
    "pbex_dalton",
    "rpbex_dalton",
    "pbex",
    "revpbex_untested",
    "pw92c_untested",
    "tfk_untested",
    "pw91k_untested",
  };

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
    case XC_R_ZI_ZJ_ZK:
      return 4;
    case XC_ABFGH:
    case XC_RSZXY:
      return 5;
    default:
      assert(0 && "Invalid xc_mode");
      return -1;
    }
}

// This routine can be used independently of the number
// of variables, i.e. LDA/GGA spin polarized/unpolarized
template<class num>
static void sum_functionals(num &res,
			    const densvars<num> &d)
{
  res = 0;

  if (xc_params[XC_VWN5_CORRELATION] != 0)
    res += xc_params[XC_VWN5_CORRELATION]*
      vwn5_correlation(d.R,d.S,d.zeta);

  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*
      slater_exchange(d.rhoa,d.rhob); 

  if (xc_params[XC_LYP_CORRELATION] != 0)
    res += xc_params[XC_LYP_CORRELATION]*
      lyp_correlation(d.R,d.rhoa,d.rhob,d.Z,d.grada2,d.gradb2);

  if (xc_params[XC_PW91_EXCHANGE] != 0)
    res += xc_params[XC_PW91_EXCHANGE]*
      pw91_exchange(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_PBE_CORRELATION] != 0)
    res += xc_params[XC_PBE_CORRELATION]*
      pbe_correlation(d.R,d.S,d.Z,d.rhoa,d.rhob,d.zeta);

 if (xc_params[XC_OLD_PBE_CORRELATION] != 0)
    res += xc_params[XC_OLD_PBE_CORRELATION]*
      pbe_correlation_unpolarized(d.R,d.Z);

  if (xc_params[XC_BECKE_EXCHANGE] != 0)
    res += xc_params[XC_BECKE_EXCHANGE]*
      becke_exchange(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_OLD_PBE_EXCHANGE] != 0)
    res += xc_params[XC_OLD_PBE_EXCHANGE]*
      pbe_exchange_dalton(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_OLD_RPBE_EXCHANGE] != 0)
    res += xc_params[XC_OLD_RPBE_EXCHANGE]*
      rpbe_exchange_dalton(d.rhoa,d.rhob,d.grada2,d.gradab,d.gradb2);

  if (xc_params[XC_PBE_EXCHANGE] != 0)
    res += xc_params[XC_PBE_EXCHANGE]*
      pbe_exchange(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_REVPBE_EXCHANGE_UNTESTED] != 0)
    res += xc_params[XC_REVPBE_EXCHANGE_UNTESTED]*
      revpbe_exchange(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_PW92_CORRELATION] != 0)
    res += xc_params[XC_PW92_CORRELATION]*
      pw92_correlation(d.R,d.zeta);

  if (xc_params[KIN_TF_UNTESTED] != 0)
    res += xc_params[KIN_TF_UNTESTED]*
      thomasfermi_kinetic(d.rhoa,d.rhob);

  if (xc_params[KIN_PW91_UNTESTED] != 0)
    res += xc_params[KIN_PW91_UNTESTED]*
      pw91_kinetic(d.rhoa,d.rhob,d.grada2,d.gradb2);

  res.deriv_facs();
}
			    
template<class num, int N>
static void setdens_rszxy(densvars<taylor<num,5,N> >&d, const num dd[5])
{
  d.R = taylor<num,5,N>(dd[0],0);
  d.S = taylor<num,5,N>(dd[1],1);
  d.Z = taylor<num,5,N>(dd[2],2);
  d.X = taylor<num,5,N>(dd[3],3);
  d.Y = taylor<num,5,N>(dd[4],4);
  d.rhoa = 0.5*(d.R+d.S);
  d.rhob = 0.5*(d.R-d.S);
  d.grada2 = 0.25*(d.Z+d.X+2*d.Y);
  d.gradab = 0.25*(d.Z-d.X);
  d.gradb2 = 0.25*(d.Z+d.X-2*d.Y);
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void setdens_abfgh(densvars<taylor<num,5,N> >&d, const num dd[5])
{
  d.rhoa = taylor<num,5,N>(dd[0],0);
  d.rhob = taylor<num,5,N>(dd[1],1);
  d.grada2 = taylor<num,5,N>(dd[2],2);
  d.gradb2 = taylor<num,5,N>(dd[3],3);
  d.gradab = taylor<num,5,N>(dd[4],4);
  d.R = d.rhoa + d.rhob;
  d.S = d.rhoa - d.rhob;
  d.Z = d.grada2+2*d.gradab+d.gradb2; 
  d.X = d.grada2-2*d.gradab+d.gradb2; 
  d.Y = d.grada2-d.gradb2;
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void setdens_rz(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.R = taylor<num,2,N>(dd[0],0);
  d.S = 0;
  d.Z = taylor<num,2,N>(dd[1],1);
  d.X = 0;
  d.Y = 0;
  d.rhoa = 0.5*d.R;
  d.rhob = d.rhoa;
  d.grada2 = 0.25*d.Z;
  d.gradab = d.grada2;
  d.gradb2 = d.grada2;
  d.zeta = 0;
}

template<class num, int Nvar, int N>
static void setdens_rz(densvars<taylor<num,Nvar,N> >&d, const taylor<num,Nvar,N> dd[2])
{
  d.R = dd[0];
  d.S = 0;
  d.Z = dd[1];
  d.X = 0;
  d.Y = 0;
  d.rhoa = 0.5*d.R;
  d.rhob = d.rhoa;
  d.grada2 = 0.25*d.Z;
  d.gradab = d.grada2;
  d.gradb2 = d.grada2;
  d.zeta = 0;
}

template<class num, int N>
static void setdens_af(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.rhoa = taylor<num,2,N>(dd[0],0);
  d.rhob = 0;
  d.grada2 = taylor<num,2,N>(dd[1],1);
  d.gradb2 = 0;
  d.gradab = 0;
  d.R = d.rhoa;
  d.S = d.rhoa;
  d.Z = d.grada2; 
  d.X = d.grada2; 
  d.Y = d.grada2;
  d.zeta = 1;
}

template<class num, int N>
static void setdens_rs(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.R = taylor<num,2,N>(dd[0],0);
  d.S = taylor<num,2,N>(dd[1],1);
  d.rhoa = 0.5*(d.R+d.S);
  d.rhob = 0.5*(d.R-d.S);
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void setdens_ab(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.rhoa = taylor<num,2,N>(dd[0],0);
  d.rhob = taylor<num,2,N>(dd[1],1);
  d.R = d.rhoa + d.rhob;
  d.S = d.rhoa - d.rhob;
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void setdens_a(densvars<taylor<num,1,N> >&d, const num dd[1])
{
  d.rhoa = taylor<num,1,N>(dd[0],0);
  d.rhob = 0;
  d.R = d.rhoa;
  d.S = d.rhoa;
  d.zeta = 1;
}

template<class num, int N>
static void setdens_r(densvars<taylor<num,1,N> >&d, const num dd[1])
{
  d.R = taylor<num,1,N>(dd[0],0);
  d.S = 0;
  d.rhoa = 0.5*d.R;
  d.rhob = 0.5*d.R;
  d.zeta = 0;
}

template<class num, int Nvar, int N>
static void setdens_r(densvars<taylor<num,Nvar,N> >&d, 
		      const taylor<num,Nvar,N> dd[1])
{
  d.R = dd[0];
  d.S = 0;
  d.rhoa = 0.5*d.R;
  d.rhob = 0.5*d.R;
  d.zeta = 0;
}

template<class num, int N>
static void setdens_r_ri_rj_rk(densvars<taylor<num,4,N> >&d, const num dd[4])
{
  d.R = taylor<num,4,N>(dd[0],0);
  d.S = 0;
  taylor<num,4,N> zi(dd[1],1),zj(dd[2],2),zk(dd[3],3);
  d.Z = zi*zi + zj*zj + zk*zk;
  d.X = 0;
  d.Y = 0;
  d.rhoa = 0.5*(d.R+d.S);
  d.rhob = 0.5*(d.R-d.S);
  d.grada2 = 0.25*(d.Z+d.X+2*d.Y);
  d.gradab = 0.25*(d.Z-d.X);
  d.gradb2 = 0.25*(d.Z+d.X-2*d.Y);
  d.zeta = d.S/d.R;
}


template<class num, int Nvar, int N>
static void setdens_r_ri_rj_rk(densvars<taylor<num,Nvar,N> >&d, const taylor<num,Nvar,N> dd[4])
{
  d.R = dd[0];
  d.S = 0;
  d.Z = dd[1]*dd[1] + dd[2]*dd[2] + dd[3]*dd[3];
  d.X = 0;
  d.Y = 0;
  d.rhoa = 0.5*(d.R+d.S);
  d.rhob = 0.5*(d.R-d.S);
  d.grada2 = 0.25*(d.Z+d.X+2*d.Y);
  d.gradab = 0.25*(d.Z-d.X);
  d.gradb2 = 0.25*(d.Z+d.X-2*d.Y);
  d.zeta = d.S/d.R;
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
      else
	setdens_a(dv,d);
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
	case XC_AF: setdens_af(dv,d); break;
	default:
	  break;
	}
      sum_functionals(res,dv);
    }
  else if (fun_nvars(current_mode) == 4)
    {
      taylor<num,4,N> &res = *reinterpret_cast<taylor<num,4,N> *>(result);
      densvars<taylor<num,4,N> > dv;
      setdens_r_ri_rj_rk(dv,d);
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
    }
  else
    {
      assert(0 && "BUG: Invalid functional type");
    }
}

void xc_eval(double *result, int order, const double densvars[])
{
  assert(order>=0);
  assert(order<=XC_MAX_ORDER);
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

template<>
void setup_order<-1>(void) {}


void setup()
{
  {
    // Test if taylor structures have padding. If they do then the
    // contraction code won't work, because it relies on that the
    // 1,1,1,1.. element of the tensor has a predictable address.
    tensored_taylor<3,double,1,1>::type t;
    assert(&t[1][1][1] - &t[0][0][0] == 7 && "taylor objects have padding, nesting them will be confusing.");
  }
  setup_order<XC_MAX_ORDER>();
  is_setup = 1;
}

#ifdef NEW_CONTRACTION

#include <iostream>
using namespace std;

// Calculates a taylor expansion of the gradient (wrt density
// variable) element grad_var, in the Npt variables given by d.
// The idea is to call this with lower order perturbed densities
// as functions of perturbation strengths (Npt perturbations). 
// The highest order coefficients of the output are then the
// right hand side of a response equation.
template<class num, int Npt, int Nvar>
struct contractor
{
  typedef typename tensored_taylor<Npt,num,1,1>::type vartype;
  typedef taylor<typename tensored_taylor<Npt,num,1,1>::type,Nvar,1> gradtype;

  static void contract(num res[Nvar], 
		       const vartype d[Nvar])
  {
    // Use the "trick" that 
    // f(a+x) = f(a) + f'(a)x
    // We want to taylor expand f'(a), so we can just pick out the term
    // linear in x. We avoid calculating higher derivatives of x by
    // defining a linear polynomial in x with polynomial coefficients a
    // and 1.  This actually gives us what we want, from the magic of
    // taylor expansions. Note that we cannot insert a as a polynomial
    // after expanding in x, we have to insert both at the same time.
    // (cf. differentiating after putting in the variable value).
    gradtype seed[Nvar];

    for (int i=0;i<Nvar;i++)
      {
	seed[i][0] = d[i];
	for (int j=0;j<Nvar;j++)
	  if (j == i)
	    seed[i][j+1] = 1;
	  else
	    seed[i][j+1] = 0;
      }


    densvars<gradtype> dv;
    setdens_rz(dv,seed);

    gradtype out;
    sum_functionals(out, dv);
    //Return the 1,1,1,1 .. element of the out tensor
    //This relies on no padding, have to check if this is true
    for (int i=0;i<Nvar;i++)
      res[i] = *((1 << Npt) - 1 + (num *)&out[i+1]);

  }
};


void test_contractions_rz(void)
{
  contractor<double,2,2>::vartype d[2];
  double res[2];

  d[0][0][0] = 0.39823175667066241;
  d[0][1][0] = 0.46570334454063117;
  d[0][0][1] = 0.46570334454063117;
  d[0][1][1] = 0;

  d[1][0][0] = 0.14893941181678338;
  d[1][1][0] = -0.17835326257152995;
  d[1][0][1] = -0.17835326257152995;
  d[1][1][1] = 0.10678806018517961;

  contractor<double,2,2>::contract(res,d);

  cout << "results:" << res[0] << " " << res[1] << endl;

}



int main(void)
{
  cout.precision(15);
  double parm[XC_PARAMS_LEN] = {0};
  parm[XC_LYP_CORRELATION] = 1;
  parm[XC_SLATER_EXCHANGE] = 1;
  parm[XC_BECKE_EXCHANGE] = 1;
  xc_set_functional(XC_RZ,parm);
  //  test_contractions_cover();
  test_contractions_rz();
  return 0;
}

#endif
