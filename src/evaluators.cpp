#include "xcfun_internal.h"
#include <cstdlib>
#include <cstdio>

// Used for regularizing input
#define TINYDENS 1e-14

// This is a replacement for copysign from C99
static inline double cpsign(double x, double y)
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

static struct evaluator_table
{
  void construct(void)
  {
    for (int i=0;i<XC_NR_MODES;i++)
      for (int j=0;j<XC_NR_TYPES;j++)
	for (int k=0;k<=XC_MAX_ORDER;k++)
	  tab[i][j][k] = 0;
  }
  evaluator tab[XC_NR_MODES][XC_NR_TYPES][XC_MAX_ORDER+1];
} *eval_tab = 0;

template<class T, class scalar>
static void sum_functionals(const double *weights,
		     //		     const array<functional *> &funs,
		     scalar *res,
		     const densvars<T> &dv)
{
  T &r = *reinterpret_cast<T *>(res);
  r = 0;
  for (int i=0;i<XC_NR_PARAMS;i++)
    if (weights[i] != 0)
      {
	functional *f = xcint_functional(i);
	if (f)
	  r += weights[i]*f->eval(dv);
      }
  r.deriv_facs();
}

// 100% spin polarized LDA
template<class T, int Ndeg>
static void eval_lda_a(const xc_functional_data &fun, 
		       T *res, const T *d)
{
  typedef taylor<T,1,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.a = ttype(d[0],0);
#else
  dv.a = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
#endif
  dv.b = 0;
  dv.n = dv.a;
  dv.s = dv.a;
  dv.zeta = 1;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = 0;
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

// unpolarized LDA
template<class T, int Ndeg>
static void eval_lda_n(const xc_functional_data &fun, 
		       T *res, const T *d)
{
  typedef taylor<T,1,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.n = ttype(d[0],0);
#else
  dv.n = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
#endif
  dv.a = dv.n/2;
  dv.b = dv.a;
  dv.s = 0;
  dv.zeta = 0;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_lda_ab(const xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,2,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.a = ttype(d[0],0);
  dv.b = ttype(d[1],1);
#else
  dv.a = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  dv.b = ttype(d[1] > TINYDENS ? d[1] : TINYDENS,1);
#endif
  dv.n = dv.a+dv.b;
  dv.s = dv.a-dv.b;
  dv.zeta = dv.s/dv.n;
  //  printf("zeta = %.15f\n",dv.zeta[0]);
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  //printf("r_s = %.15f\n",dv.r_s[0]);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_lda_ns(const xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,2,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.n = ttype(d[0],0);
  dv.s = ttype(d[1],1);
#else
  dv.n = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  // Not differentiable at complete spin polarization
  dv.s = ttype(fabs(d[1]) < d[0] ? d[1] : cpsign(d[0]-TINYDENS,d[1]),1);
#endif
  dv.a = 0.5*(dv.n + dv.s);
  dv.b = 0.5*(dv.n - dv.s);
  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_gga_ab(const xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,5,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.a = ttype(d[0],0);
  dv.b = ttype(d[1],1);
#else
  dv.a = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  dv.b = ttype(d[1] > TINYDENS ? d[1] : TINYDENS,1);
#endif  

  dv.n = dv.a+dv.b;
  dv.s = dv.a-dv.b;
    
#ifdef XC_NO_REGULARIZATION
  dv.gaa = ttype(d[2],2);
  dv.gab = ttype(d[3],3);
  dv.gbb = ttype(d[4],4);
#else
  dv.gaa = ttype(d[2] >= 0 ? d[2] : 0,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
  dv.gab = ttype(d[3],3); 
#else
  dv.gab = ttype(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
  dv.gbb = ttype(d[4] >= 0 ? d[4] : 0,4);
#endif
  
  dv.gnn  = dv.gaa + 2*dv.gab + dv.gbb; 
  dv.gss  = dv.gaa - 2*dv.gab + dv.gbb;
  dv.gns  = dv.gaa - dv.gbb;
  
  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_gga_ns(const xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,5,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.n = ttype(d[0],0);
  dv.s = ttype(d[1],1);
#else
  dv.n = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  // Not differentiable at complete spin polarization
  dv.s = ttype(fabs(d[1]) < d[0] ? d[1] : cpsign(d[0]-TINYDENS,d[1]),1);
#endif
  dv.a = (dv.n + dv.s)/2;
  dv.b = (dv.n - dv.s)/2;
    
#ifdef XC_NO_REGULARIZATION
  dv.gnn = ttype(d[2],2);
  dv.gns = ttype(d[3],3);
  dv.gss = ttype(d[4],4);
#else
  dv.gnn = ttype(d[2] >= 0 ? d[2] : 0,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
  dv.gns = ttype(d[3], 3); 
#else
  dv.gns = ttype(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
  dv.gss = ttype(d[4] >= 0 ? d[4] : 0,4);
#endif  

  dv.gaa  = (dv.gss+2*dv.gns+dv.gnn)/4;
  dv.gbb  = (dv.gss-2*dv.gns+dv.gnn)/4;
  dv.gab  = (dv.gnn-dv.gss)/4;

  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}


template<class T, int Ndeg>
static void eval_gga_n(const xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,2,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.n = ttype(d[0],0);
#else
  dv.n = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
#endif
  dv.s = 0;
  dv.a = (dv.n + dv.s)/2;
  dv.b = (dv.n - dv.s)/2;
    
#ifdef XC_NO_REGULARIZATION
  dv.gnn = ttype(d[1],1);
#else
  dv.gnn = ttype(d[1] >= 0 ? d[1] : 0, 1);
#endif  
  dv.gns = 0;
  dv.gss = 0;

  dv.gaa  = (dv.gss+2*dv.gns+dv.gnn)/4;
  dv.gbb  = (dv.gss-2*dv.gns+dv.gnn)/4;
  dv.gab  = (dv.gnn-dv.gss)/4;

  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_mgga_ab(const xc_functional_data &fun, 
			 T *res, const T *d)
{
  typedef taylor<T,7,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.a = ttype(d[0],0);
  dv.b = ttype(d[1],1);
#else
  dv.a = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  dv.b = ttype(d[1] > TINYDENS ? d[1] : TINYDENS,1);
#endif  

  dv.n = dv.a+dv.b;
  dv.s = dv.a-dv.b;
    
#ifdef XC_NO_REGULARIZATION
  dv.gaa = ttype(d[2],2);
  dv.gab = ttype(d[3],3);
  dv.gbb = ttype(d[4],4);
#else
  dv.gaa = ttype(d[2] >= 0 ? d[2] : 0,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
  dv.gab = ttype(d[3],3); 
#else
  dv.gab = ttype(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
  dv.gbb = ttype(d[4] >= 0 ? d[4] : 0,4);
#endif  

  dv.gnn  = dv.gaa + 2*dv.gab + dv.gbb; 
  dv.gss  = dv.gaa - 2*dv.gab + dv.gbb;
  dv.gns  = dv.gaa - dv.gbb;
  
#ifdef XC_NO_REGULARIZATION
  dv.taua = ttype(d[5], 5);
  dv.taub = ttype(d[6], 6);
#else
  dv.taua = ttype(d[5] >= TINYDENS ? d[5] : TINYDENS, 5);
  dv.taub = ttype(d[6] >= TINYDENS ? d[6] : TINYDENS, 6);
#endif
  dv.tau = dv.taua + dv.taub;

  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_mgga_ns(const xc_functional_data &fun, 
			 T *res, const T *d)
{
  typedef taylor<T,7,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.n = ttype(d[0],0);
  dv.s = ttype(d[1],1);
#else
  dv.n = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  // Not differentiable at complete spin polarization
  dv.s = ttype(fabs(d[1]) < d[0] ? d[1] : cpsign(d[0]-TINYDENS,d[1]),1);
#endif
  dv.a = (dv.n + dv.s)/2;
  dv.b = (dv.n - dv.s)/2;
    
#ifdef XC_NO_REGULARIZATION
  dv.gnn = ttype(d[2],2);
  dv.gns = ttype(d[3],3);
  dv.gss = ttype(d[4],4);
#else
  dv.gnn = ttype(d[2] >= 0 ? d[2] : TINYDENS,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
  dv.gns = ttype(d[3],3); 
#else
  dv.gns = ttype(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
  dv.gss = ttype(d[4] >= 0 ? d[4] : 0,4);
#endif
  
  dv.gaa  = (dv.gss+2*dv.gns+dv.gnn)/4;
  dv.gbb  = (dv.gss-2*dv.gns+dv.gnn)/4;
  dv.gab  = (dv.gnn-dv.gss)/4;

#ifdef XC_NO_REGULARIZATION
  dv.tau = ttype(d[5],5);
  dv.taua = 0.5*(d[5]+ttype(d[6],6));
  dv.taub = dv.tau - dv.taua;
#else
  dv.tau = ttype(d[5] >= TINYDENS ? d[5] : TINYDENS, 5);
  dv.taua = 0.5*(dv.tau + 
		 ttype(fabs(d[6]) < dv.tau ? 
		       d[6] : cpsign(dv.tau[0]-TINYDENS,d[6]), 6));
  dv.taub = dv.tau - dv.taua;
#endif

  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_mlgga_ab(const xc_functional_data &fun, 
			 T *res, const T *d)
{
  typedef taylor<T,9,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.a = ttype(d[0],0);
  dv.b = ttype(d[1],1);
#else
  dv.a = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  dv.b = ttype(d[1] > TINYDENS ? d[1] : TINYDENS,1);
#endif  

  dv.n = dv.a+dv.b;
  dv.s = dv.a-dv.b;
    
#ifdef XC_NO_REGULARIZATION
  dv.gaa = ttype(d[2],2);
  dv.gab = ttype(d[3],3);
  dv.gbb = ttype(d[4],4);
#else
  dv.gaa = ttype(d[2] >= 0 ? d[2] : 0,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
  dv.gab = ttype(d[3],3); 
#else
  dv.gab = ttype(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
  dv.gbb = ttype(d[4] >= 0 ? d[4] : 0,4);
#endif  

  dv.gnn  = dv.gaa + 2*dv.gab + dv.gbb; 
  dv.gss  = dv.gaa - 2*dv.gab + dv.gbb;
  dv.gns  = dv.gaa - dv.gbb;
  
#ifdef XC_NO_REGULARIZATION
  dv.taua = ttype(d[5], 5);
  dv.taub = ttype(d[6], 6);
#else
  dv.taua = ttype(d[5] >= TINYDENS ? d[5] : TINYDENS, 5);
  dv.taub = ttype(d[6] >= TINYDENS ? d[6] : TINYDENS, 6);
#endif
  dv.tau = dv.taua + dv.taub;

  dv.lapa = ttype(d[7],7);
  dv.lapb = ttype(d[8],8);

  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

template<class T, int Ndeg>
static void eval_mlgga_ns(const xc_functional_data &fun, 
			  T *res, const T *d)
{
  typedef taylor<T,9,Ndeg> ttype;
  densvars<ttype> dv(fun.parameters);
#ifdef XC_NO_REGULARIZATION
  dv.n = ttype(d[0],0);
  dv.s = ttype(d[1],1);
#else
  dv.n = ttype(d[0] > TINYDENS ? d[0] : TINYDENS,0);
  // Not differentiable at complete spin polarization
  dv.s = ttype(fabs(d[1]) < d[0] ? d[1] : cpsign(d[0]-TINYDENS,d[1]),1);
#endif
  dv.a = (dv.n + dv.s)/2;
  dv.b = (dv.n - dv.s)/2;
    
#ifdef XC_NO_REGULARIZATION
  dv.gnn = ttype(d[2],2);
  dv.gns = ttype(d[3],3);
  dv.gss = ttype(d[4],4);
#else
  dv.gnn = ttype(d[2] >= 0 ? d[2] : TINYDENS,2);
#ifdef XC_NO_SCHWARZ_REGULARIZATION
  dv.gns = ttype(d[3],3); 
#else
  dv.gns = ttype(d[3]*d[3] <= d[2]*d[4] ? d[3] : sqrt(d[2]*d[4]),3); 
#endif
  dv.gss = ttype(d[4] >= 0 ? d[4] : 0,4);
#endif
  
  dv.gaa  = (dv.gss+2*dv.gns+dv.gnn)/4;
  dv.gbb  = (dv.gss-2*dv.gns+dv.gnn)/4;
  dv.gab  = (dv.gnn-dv.gss)/4;

#ifdef XC_NO_REGULARIZATION
  dv.tau = ttype(d[5],5);
  dv.taua = 0.5*(d[5]+ttype(d[6],6));
  dv.taub = dv.tau - dv.taua;
#else
  dv.tau = ttype(d[5] >= TINYDENS ? d[5] : TINYDENS, 5);
  dv.taua = 0.5*(dv.tau + 
		 ttype(fabs(d[6]) < dv.tau ? 
		       d[6] : cpsign(dv.tau[0]-TINYDENS,d[6]), 6));
  dv.taub = dv.tau - dv.taua;
#endif

  dv.lapa = ttype(0.5*(d[7]+d[8]),7);
  dv.lapb = ttype(-0.5*(d[7]-d[8]),8);

  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  dv.n_m13 = pow(dv.n,-1.0/3.0);
  dv.a_43 = pow(dv.a,4.0/3.0);
  dv.b_43 = pow(dv.b,4.0/3.0);
  sum_functionals(fun.parameters,/* fun.active_functionals, */res,dv);
}

// Template loops to set up evaluators for all classes

template<int Ndeg>
static void eval_setup_lda(void)
{
  eval_tab->tab[XC_VARS_A][XC_LDA][Ndeg] = eval_lda_a<double,Ndeg>;
  eval_tab->tab[XC_VARS_AB][XC_LDA][Ndeg] = eval_lda_ab<double,Ndeg>;
  eval_tab->tab[XC_VARS_N][XC_LDA][Ndeg] = eval_lda_n<double,Ndeg>;
  eval_tab->tab[XC_VARS_NS][XC_LDA][Ndeg] = eval_lda_ns<double,Ndeg>;
  eval_setup_lda<Ndeg-1>();
}
template<> void eval_setup_lda<-1>(void) {}

template<int Ndeg>
static void eval_setup_gga(void)
{
  eval_tab->tab[XC_VARS_AB][XC_GGA][Ndeg] = eval_gga_ab<double,Ndeg>;
  eval_tab->tab[XC_VARS_NS][XC_GGA][Ndeg] = eval_gga_ns<double,Ndeg>;
  eval_tab->tab[XC_VARS_N][XC_GGA][Ndeg] = eval_gga_n<double,Ndeg>;
  eval_setup_gga<Ndeg-1>();
}
template<> void eval_setup_gga<-1>(void) {}

template<int Ndeg>
static void eval_setup_mgga(void)
{
  eval_tab->tab[XC_VARS_AB][XC_MGGA][Ndeg] = eval_mgga_ab<double,Ndeg>;
  eval_tab->tab[XC_VARS_NS][XC_MGGA][Ndeg] = eval_mgga_ns<double,Ndeg>;
  eval_setup_mgga<Ndeg-1>();
}
template<> void eval_setup_mgga<-1>(void) {}

template<int Ndeg>
static void eval_setup_mlgga(void)
{
  eval_tab->tab[XC_VARS_AB][XC_MLGGA][Ndeg] = eval_mlgga_ab<double,Ndeg>;
  eval_tab->tab[XC_VARS_NS][XC_MLGGA][Ndeg] = eval_mlgga_ns<double,Ndeg>;
  eval_setup_mlgga<Ndeg-1>();
}
template<> void eval_setup_mlgga<-1>(void) {}


evaluator xc_evaluator_lookup(int mode, int type, int order)
{
  if (!eval_tab)
    {
      eval_tab = (evaluator_table *)malloc(sizeof(evaluator_table));
      eval_tab->construct();
      eval_setup_lda<XC_LDA_MAX_ORDER>();
      eval_setup_gga<XC_GGA_MAX_ORDER>();
      eval_setup_mgga<XC_MGGA_MAX_ORDER>();
      eval_setup_mlgga<XC_MLGGA_MAX_ORDER>();
    }
  assert(mode>=0 && mode < XC_NR_MODES);
  assert(type>=0 && type < XC_NR_TYPES);
  assert(order>=0 && order <= XC_MAX_ORDER);
  return eval_tab->tab[mode][type][order];
}
