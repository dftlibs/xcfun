#ifndef FUNCTIONAL_H
#define FUNCTIONAL_H

// Everything required to implement a functional
#include "taylor.h"
#include "xcfun.h"
#include "config.h"
#include "array.h"
#include "specmath.h"
#include "parameters.h"

// MGGA's have 7 variables (GGA+tau)
#define XC_MAX_NVAR 7

// Use this type to make it clear what is a hard-coded parameter
// (for example for grepping for parameters in functionals)

typedef double parameter;


// Functions for _user definable parameters_, for example the range
// separation parameter mu.
//setting xc_param_lookup(const xc_functional_data *params, 
//			const char *name);
//double xc_param_get(const xc_functional::xc_functional_data *params, 
//		    const setting &s);


// Variables for expressing functionals, these are redundant because
// different functionals have different needs.
template<class T>
struct densvars
{
  densvars(const parameter *p) : params(p) {}
  //For getting user defined parameters
  const parameter *params;
  double get(enum xc_parameters p) const
  {
    assert(p>=0);
    assert(p<XC_NR_PARAMS);
    return params[p];
  }
  
  T a, b, gaa, gab, gbb;
  /* na+nb, na-nb, (grad n)^2, (grad n).(grad s), (grad s)^2 */
  T n, s, gnn, gns, gss; 

  T tau, taua, taub; // Kinetic energy densities.

  T zeta; //s/n
  T r_s; // (3/4pi)^1/3*n^(-1/3)
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
  void energy_fun(taylor<double,Nvar,Ndeg> 
		  (*f)(const densvars< taylor<double,Nvar,Ndeg> > &))
  {
    assert(Nvar<=XC_MAX_NVAR);
    assert(Ndeg<=XC_MAX_ORDER);
    ftab[Nvar][Ndeg] = reinterpret_cast<void *>(f);
  }
  template<int Nvar, int Ndeg>
  taylor<double,Nvar,Ndeg> 
  eval(const densvars< taylor<double,Nvar,Ndeg> > &dv)
  {
    assert(Nvar<=XC_MAX_NVAR);
    assert(Ndeg<=XC_MAX_ORDER);
    assert(ftab[Nvar][Ndeg]);
    taylor<double,Nvar,Ndeg> 
      (*f)(const densvars< taylor<double,Nvar,Ndeg> > &) = 
      reinterpret_cast<taylor<double,Nvar,Ndeg> 
      (*)(const densvars< taylor<double,Nvar,Ndeg> > &)>(ftab[Nvar][Ndeg]);
    return f(dv);
  }
  //protected:
  enum xc_parameters m_name;

  int m_type;
  void *ftab[XC_MAX_NVAR+1][XC_MAX_ORDER+1];

  double *test_input;
  double *test_output;
  int test_mode;
  int test_order;
  double test_threshold;
};

array<functional *> &xc_get_functional_array(void);
functional *xc_get_functional_by_name(const char *name);

// Create a functional object and run f to fill it with data
int xc_run_functional_setup(void (*f)(functional &));

//  Macros to help set up functionals to the order defined in config.h

#define SET_LDA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<1,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<2,ORDER>(ENERGY);

#define SET_GGA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<2,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<5,ORDER>(ENERGY); 

#define SET_MGGA_N(FUNOBJ,ENERGY,ORDER)\
  FUNOBJ.energy_fun<3,ORDER>(ENERGY);\
  FUNOBJ.energy_fun<7,ORDER>(ENERGY);

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

#define SET_LDA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_LDA); SETN(XC_LDA_MAX_ORDER,LDA,FUNOBJ,FUNCTION)
#define SET_GGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_GGA); SETN(XC_GGA_MAX_ORDER,GGA,FUNOBJ,FUNCTION)
#define SET_MGGA_ENERGY_FUNCTION(FUNOBJ,FUNCTION) assert(FUNOBJ.m_type <= XC_MGGA); SETN(XC_MGGA_MAX_ORDER,MGGA,FUNOBJ,FUNCTION)

#endif
