#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

// Everything required to implement a functional
#include "taylor.h"
#ifdef XCFUN_CONTRACTIONS
#include "ctaylor.h"
#endif
#define XCFUN_INTERNAL
#include "xcfun.h"
#include "config.h"
#include "specmath.h"
#include "parameters.h"
#include <stdio.h>

// MLGGA's have 9 variables (GGA+tau+laplacian)
#define XC_MAX_NVAR 9

// Use this type to make it clear what is a hard-coded parameter
// (for example for grepping for parameters in functionals)

typedef ireal_t parameter;

// Variables for expressing functionals, these are redundant because
// different functionals have different needs.
template<class T>
struct densvars
{
  densvars(const double *p) : params(p) {}
  //For getting user defined parameters
  const double *params;
  double get_param(int param) const
  {
    assert(param>=0);
    assert(param<XC_NR_PARAMS);
    return params[param];
  }
  
  T a, b, gaa, gab, gbb;
  /* na+nb, na-nb, (grad n)^2, (grad n).(grad s), (grad s)^2 */
  T n, s, gnn, gns, gss; 

  T tau, taua, taub; // Kinetic energy densities.

  T lapa,lapb; // Density Laplacians

  T zeta; //s/n
  T r_s; // (3/4pi)^1/3*n^(-1/3)
  T n_m13; // pow(n,-1.0/3.0)
  T a_43, b_43; // pow(a,4.0/3.0), pow(b,4.0/3.0)
};

class functional
{
public:
  void construct();
  void destroy();
  void describe(enum xc_parameters weight_param, int type, const char *oneliner,
		const char *reference);
  void describe_parameter(enum xc_parameters param, const char *description,
			  double default_value);
  void add_test(int mode, int order,
		const double *test_in,
		const double *test_out,
		double threshold);
  int validate();
  template<int Nvar, int Ndeg>
  void energy_fun(taylor<ireal_t,Nvar,Ndeg> 
		  (*f)(const densvars< taylor<ireal_t,Nvar,Ndeg> > &))
  {
    assert(Nvar<=XC_MAX_NVAR);
    assert(Ndeg<=XC_MAX_ORDER);
    // See comment at eval about this cast
#ifdef _MSC_VER
    ftab[Nvar][Ndeg] = *reinterpret_cast<void **>(&f);
#else
    ftab[Nvar][Ndeg] = reinterpret_cast<void *>(f);
#endif
  }
#ifdef XCFUN__CONTRACTIONS
  template<int Ndeg>
  void contraction_fun(ctaylor<ireal_t,Ndeg> 
		       (*f)(const densvars< ctaylor<ireal_t,Ndeg> > &))
  {
    assert(Ndeg<=XC_CONTRACT_MAX_ORDER);
    contract_ftab[Ndeg] = reinterpret_cast<void *>(f);
  }
#endif
  template<int Nvar, int Ndeg>
  taylor<ireal_t,Nvar,Ndeg> 
  eval(const densvars< taylor<ireal_t,Nvar,Ndeg> > &dv)
  {
    assert(Nvar<=XC_MAX_NVAR);
    assert(Ndeg<=XC_MAX_ORDER);
    assert(ftab[Nvar][Ndeg]);
#ifdef _MSC_VER
    /* Here we need to cast the void pointer back to a function
       pointer of the appropriate type. This is supposed to work
       under C++0x, may not be supported under all compilers.
       Using solution from 
       http://stackoverflow.com/questions/1096341/function-pointers-casting-in-c
     */
    taylor<ireal_t,Nvar,Ndeg> 
      (*f)(const densvars< taylor<ireal_t,Nvar,Ndeg> > &);
    *reinterpret_cast<void**>(&f) = ftab[Nvar][Ndeg];
#else
    taylor<ireal_t,Nvar,Ndeg> 
      (*f)(const densvars< taylor<ireal_t,Nvar,Ndeg> > &) = 
      reinterpret_cast<taylor<ireal_t,Nvar,Ndeg> 
      (*)(const densvars< taylor<ireal_t,Nvar,Ndeg> > &)>(ftab[Nvar][Ndeg]);
#endif
    return f(dv);
  }
#ifdef XCFUN_CONTRACTIONS
  template<int Ndeg>  
  ireal_t contract_eval(const densvars< ctaylor<ireal_t,Ndeg> > &dv)
  {
    assert(Ndeg<=XC_CONTRACT_MAX_ORDER);
    assert(contract_ftab[Ndeg]);
    ctaylor<ireal_t,Ndeg> 
      (*f)(const densvars< ctaylor<ireal_t,Ndeg> > &) = 
      reinterpret_cast<ctaylor<ireal_t,Ndeg> 
      (*)(const densvars< ctaylor<ireal_t,Ndeg> > &)>(contract_ftab[Ndeg]);
    return f(dv);
  }
#endif

  //protected:
  enum xc_parameters m_name;

  int m_type;
  void *ftab[XC_MAX_NVAR+1][XC_MAX_ORDER+1];
#ifdef XCFUN_CONTRACTIONS
  void *contract_ftab[XC_CONTRACT_MAX_ORDER+1]; 
#endif
  double *test_input;
  double *test_output;
  int test_mode;
  int test_order;
  double test_threshold;
};

functional *xc_get_functional_by_name(const char *name);

// Create a functional object and run f to fill it with data
int xc_run_functional_setup(void (*f)(functional &));

//  Macros to help set up functionals to the order defined in config.h

#ifdef XCFUN_CONTRACTIONS
#define SET_CONTRACT_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.contraction_fun<ORDER>(ENERGY);
#endif

#define SET_LDA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<1,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<2,ORDER>(ENERGY);

#define SET_GGA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<2,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<5,ORDER>(ENERGY); 

#define SET_MGGA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<3,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<7,ORDER>(ENERGY);

#define SET_MLGGA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<4,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<9,ORDER>(ENERGY);


//Extend this list if you want higher order than 9
#define SETN0(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,0)
#define SETN1(TYPE,F,ERG) SETN0(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,1)
#define SETN2(TYPE,F,ERG) SETN1(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,2)
#define SETN3(TYPE,F,ERG) SETN2(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,3)
#define SETN4(TYPE,F,ERG) SETN3(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,4)
#define SETN5(TYPE,F,ERG) SETN4(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,5)
#define SETN6(TYPE,F,ERG) SETN5(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,6)
#define SETN7(TYPE,F,ERG) SETN6(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,7)
#define SETN8(TYPE,F,ERG) SETN7(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,8)
#define SETN9(TYPE,F,ERG) SETN8(TYPE,F,ERG) SET_##TYPE##_N(F,ERG,9)

// Since we cannot use -1 to say that we want no items at all, use 666 instead
#define SETN666(TYPE,F,ERG)

// cpp quirk to force expansion of N before pasting 
#define SETNQ(N,TYPE,F,ERG) SETN ## N (TYPE,F,ERG)
#define SETN(N,TYPE,F,ERG) SETNQ(N,TYPE,F,ERG)

#ifdef XCFUN_CONTRACTIONS
#define SET_LDA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_LDA); \
SETN(XC_LDA_MAX_ORDER,LDA,FUNOBJ,FUNCTION)\
SETN(XC_CONTRACT_MAX_ORDER,CONTRACT,FUNOBJ,FUNCTION)

#define SET_GGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_GGA); \
SETN(XC_GGA_MAX_ORDER,GGA,FUNOBJ,FUNCTION)\
SETN(XC_CONTRACT_MAX_ORDER,CONTRACT,FUNOBJ,FUNCTION)

#define SET_MGGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_MGGA); \
SETN(XC_MGGA_MAX_ORDER,MGGA,FUNOBJ,FUNCTION)\
SETN(XC_CONTRACT_MAX_ORDER,CONTRACT,FUNOBJ,FUNCTION)

#define SET_MLGGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_MLGGA); \
SETN(XC_MLGGA_MAX_ORDER,MLGGA,FUNOBJ,FUNCTION)\
SETN(XC_CONTRACT_MAX_ORDER,CONTRACT,FUNOBJ,FUNCTION)

#else
#define SET_LDA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_LDA); \
SETN(XC_LDA_MAX_ORDER,LDA,FUNOBJ,FUNCTION)

#define SET_GGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_GGA); \
SETN(XC_GGA_MAX_ORDER,GGA,FUNOBJ,FUNCTION)

#define SET_MGGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_MGGA); \
SETN(XC_MGGA_MAX_ORDER,MGGA,FUNOBJ,FUNCTION)

#define SET_MLGGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_MLGGA); \
SETN(XC_MLGGA_MAX_ORDER,MLGGA,FUNOBJ,FUNCTION)
#endif

#endif
