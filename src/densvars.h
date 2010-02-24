#ifndef DENSVARS_H
#define DENSVARS_H

#include "xc_fun.h"
#include "taylor.h"

template<class T> 
struct densvars
{
  // Variables for expressing functionals, these are redundant because
  // different functionals have different needs.

  /* na+nb, na-nb, (grad n)^2, (grad n).(grad s), (grad s)^2 */
  T n, s, gnn, gns, gss; 
  T a, b, gaa, gab, gbb;
#ifdef METAGGA
  T tau, taua, taub;
#endif
  T zeta; //s/n
  T r_s; // (3/4pi)^1/3*n^(-1/3)
};

template<int mode>
struct xc_mode_info
{
  enum { Nvar = 0 };
};

template<>
struct xc_mode_info<XC_ALDA>
{
  enum { Nvar = 1 };
};

template<>
struct xc_mode_info<XC_RLDA>
{
  enum { Nvar = 1 };
};

template<>
struct xc_mode_info<XC_ABLDA>
{
  enum { Nvar = 2 };
};

template<>
struct xc_mode_info<XC_RSLDA>
{
  enum { Nvar = 2 };
};

template<>
struct xc_mode_info<XC_AGGA>
{
  enum { Nvar = 2 };
};

template<>
struct xc_mode_info<XC_RGGA>
{
  enum { Nvar = 2 };
};

template<>
struct xc_mode_info<XC_ABGGA>
{
  enum { Nvar = 5 };
};

template<>
struct xc_mode_info<XC_RSGGA>
{
  enum { Nvar = 5 };
};


template<class num, enum xc_mode mode, int Ndeg>
static void setdens(densvars<taylor<num,xc_mode_info<mode>::Nvar,Ndeg> >&d, 
		    const num *dd)
{
  typedef taylor<num,xc_mode_info<mode>::Nvar,Ndeg> ttype;
  switch (mode)
    {
    case XC_RSGGA:
      d.n = ttype(dd[0],0);
      d.s = ttype(dd[1],1);
      d.gnn = ttype(dd[2],2);
      d.gns = ttype(dd[3],3);
      d.gss = ttype(dd[4],4);
      d.a = 0.5*(d.n+d.s);
      d.b = 0.5*(d.n-d.s);
      d.gaa = 0.25*(d.gnn+d.gns+2*d.gss);
      d.gab = 0.25*(d.gnn-d.gns);
      d.gbb = 0.25*(d.gnn+d.gns-2*d.gss);
      break;
    case XC_ABGGA:
      d.a = ttype(dd[0],0);
      d.b = ttype(dd[1],1);
      d.gaa = ttype(dd[2],2);
      d.gbb = ttype(dd[3],3);
      d.gab = ttype(dd[4],4);
      d.n = d.a + d.b;
      d.s = d.a - d.b;
      d.gnn = d.gaa+2*d.gab+d.gbb; 
      d.gns = d.gaa-2*d.gab+d.gbb;
      d.gss = d.gaa-d.gbb;
      break;
    case XC_RGGA:
      d.n = ttype(dd[0],0);
      d.s = 0;
      d.gnn = ttype(dd[1],1);
      d.gns = 0;
      d.gss = 0;
      d.a = 0.5*d.n;
      d.b = 0.5*d.n;
      d.gaa = 0.25*d.gnn;
      d.gab = 0.25*d.gnn;
      d.gbb = 0.25*d.gnn;
    default:
      assert(0 && "Unimplemented mode");
      break;
    };
  d.zeta = d.s/d.n;
  d.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(d.n,-1.0/3.0);
}
#if 0
template<class num, int N>
static void setdens_rs(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.n = taylor<num,2,N>(dd[0],0);
  d.s = taylor<num,2,N>(dd[1],1);
  d.a = 0.5*(d.n+d.s);
  d.b = 0.5*(d.n-d.s);
  d.zeta = d.s/d.n;
  d.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(d.n,-1.0/3.0);
}

template<class num, int N>
static void setdens_ab(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.a = taylor<num,2,N>(dd[0],0);
  d.b = taylor<num,2,N>(dd[1],1);
  d.n = d.a + d.b;
  d.s = d.a - d.b;
  d.zeta = d.s/d.n;
  d.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(d.n,-1.0/3.0);
}

template<class num, int N>
static void setdens_r(densvars<taylor<num,1,N> >&d, const num dd[1])
{
  d.n = taylor<num,1,N>(dd[0],0);
  d.s = 0;
  d.a = 0.5*d.n;
  d.b = 0.5*d.n;
  d.zeta = 0;
  d.r_s = pow(3/(4*M_PI),1.0/3.0)*pow(d.n,-1.0/3.0);
}
#endif

#endif
