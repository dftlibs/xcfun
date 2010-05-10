
template<class T, class scalar>
void sum_functionals(const array<double> &weights, 
		     const array<functional *> &funs,
		     scalar *res,
		     const densvars<T> &dv)
{
  T &r = *reinterpret_cast<T *>(res);
  r = 0;
  for (int i=0;i<weights.size();i++)
    r += weights[i]*funs[i]->eval(dv);
  r.deriv_facs();
}

// 100% spin polarized LDA
template<class T, int Ndeg>
static void eval_lda_a(const xc_functional::xc_functional_data &fun, 
		       T *res, const T *d)
{
  typedef taylor<T,1,Ndeg> ttype;
  densvars<ttype> dv(&fun);
  dv.a = ttype(d[0],0);
  dv.b = 0;
  dv.n = dv.a;
  dv.s = dv.a;
  dv.zeta = 1;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  sum_functionals(fun.weights,fun.active_functionals,res,dv);
}

// unpolarized LDA
template<class T, int Ndeg>
static void eval_lda_r(const xc_functional::xc_functional_data &fun, 
		       T *res, const T *d)
{
  typedef taylor<T,1,Ndeg> ttype;
  densvars<ttype> dv(&fun);
  dv.n = ttype(d[0],0);
  dv.a = dv.n/2;
  dv.b = dv.a;
  dv.s = 0;
  dv.zeta = 0;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  sum_functionals(fun.weights,fun.active_functionals,res,dv);
}

template<class T, int Ndeg>
static void eval_lda_ab(const xc_functional::xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,2,Ndeg> ttype;
  densvars<ttype> dv(&fun);
  dv.a = ttype(d[0],0);
  dv.b = ttype(d[1],1);
  dv.n = dv.a+dv.b;
  dv.s = dv.a-dv.b;
  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  sum_functionals(fun.weights,fun.active_functionals,res,dv);
}

template<class T, int Ndeg>
static void eval_gga_ab(const xc_functional::xc_functional_data &fun, 
			T *res, const T *d)
{
  typedef taylor<T,5,Ndeg> ttype;
  densvars<ttype> dv(&fun);
  dv.a = ttype(d[0],0);
  dv.b = ttype(d[1],1);
  
  dv.n = dv.a+dv.b;
  dv.s = dv.a-dv.b;
    
  dv.gaa = ttype(d[2],2);
  dv.gbb = ttype(d[3],3);
  dv.gab = ttype(d[4],4);
  
  dv.gnn  = dv.gaa + 2*dv.gab + dv.gbb; 
  dv.gss  = dv.gaa - 2*dv.gab + dv.gbb;
  dv.gns  = dv.gaa - dv.gbb;
  
  dv.zeta = dv.s/dv.n;
  dv.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(dv.n,-1.0/3.0);
  sum_functionals(fun.weights,fun.active_functionals,res,dv);
}

template<Ndeg>
void eval_setup_lda(void)
{
  
}
