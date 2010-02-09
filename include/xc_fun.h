#ifndef XC_FUN
#define XC_FUN

/* 
   Main entry point to DFT exchange-correlation functionals 
   by ulfek 2009 
*/

// Note that these parameters have be in sync with those in
// xc_fun_module.f90 (Fortran version)

// Highest order of derivatives. Note that the generated code
// becomes huge if this is too large
#define XC_MAX_ORDER 6

/*
  The following variables are used in the code:
  A = \rho_\alpha
  B = \rho_\beta
  F = \nabla A \cdot \nabla A (sometimes called \sigma_\alpha\alpha)
  G = \nabla B \cdot \nabla B (\sigma_\beta\beta)
  H = \nabla A \cdot \nabla B (\sigma_\alpha\beta)

  Alternatively you can use;
  R = \rho_\alpha + \rho_\beta
  S = \rho_\alpha - \rho_\beta
  Z = \nabla R \cdot \nabla R (sometimes called \sigma)
  X = \nabla S \cdot \nabla S
  Y = \nabla R \cdot \nabla S
*/

enum xc_mode
  {
    XC_A = 0, // Fully spin polarized LDA
    XC_R, // Unpolarized LDA
    XC_AB, // LDA with alpha/beta densities
    XC_RS, // LDA with total density/polarization density
    XC_AF, // Fully spin polarized GGA
    XC_RZ, // Unpolarized GGA
    XC_ABFGH, // Generic GGA alpha,beta variables 
    XC_RSZXY, // Generic GGA R,S variables
    XC_R_ZI_ZJ_ZK, // R and the individual gradient elements
  };

enum xc_funparams
  {
    XC_VWN5_CORRELATION = 0,
    XC_SLATER_EXCHANGE,
    XC_LYP_CORRELATION,
    XC_PW91_EXCHANGE,
    XC_PBE_CORRELATION,
    XC_OLD_PBE_CORRELATION,
    XC_BECKE_EXCHANGE,
    XC_OLD_PBE_EXCHANGE,
    XC_OLD_RPBE_EXCHANGE,
    XC_PBE_EXCHANGE,
    XC_REVPBE_EXCHANGE_UNTESTED,
    XC_PW92_CORRELATION,
    KIN_TF_UNTESTED,       // (orbital-free) kinetic energy functionals next
    KIN_PW91_UNTESTED,
    XC_PARAMS_LEN // This has to be last!
  };

int xc_parse_functional(const char *functional_name);

// Set the functional as defined by the parameters in the
// params array. 
void xc_set_functional(enum xc_mode mode, 
		       const double params[XC_PARAMS_LEN]);

//Return the total number of terms calculated by xc_eval
int xc_len(enum xc_mode mode, int order); 

// Evaluate the exchange-correlation functional defined by
// a previous call to xc_set_functional() with the density
// parameters in densvar, with derivatives up to the given
// order. The result is put in the result array, in graded
// lexicographical ordering. This function also has a fortran 
// version.
void xc_eval(double result[], int order, const double densvars[]);
#ifdef WITH_SINGLE
void xc_eval_single(float *result, int order, const float densvars[]);
#endif
//Return the index (into result) for a given derivative, i.e.
//{0,2,1,0} gives the d^3/dS^2dZ term index.
int xc_index(enum xc_mode mode, const int derivative[]);

// Return a string describing the functional and its implementation.
const char *xc_fun_reference(enum xc_funparams functions);

#ifdef WITH_QD
#include <qd/qd_real.h>
void xc_eval_qd(qd_real *result, int order, const qd_real densvars[]);
#endif

#endif
