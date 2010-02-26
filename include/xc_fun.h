#ifndef XC_FUN
#define XC_FUN

/* 
   Main entry point to DFT exchange-correlation functionals 
   by ulfek 2009 
*/

#define XC_FUN_VERSION "1.0"

// Note that these parameters have be in sync with those in
// xc_fun_module.f90 (Fortran version)

// Highest order of derivatives. Note that the generated code
// becomes huge if this is too large
#define XC_MAX_ORDER 2

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

// Which variables to use, each mode has its own bit.
enum xc_mode
  {
    XC_ALDA    = 1, // Fully spin polarized LDA
    XC_RLDA    = 2, // Unpolarized LDA
    XC_ABLDA   = 4, // LDA with alpha/beta densities
    XC_RSLDA   = 8, // LDA with total density/polarization density
    XC_AGGA    =16, // Fully spin polarized GGA
    XC_RGGA    =32, // Unpolarized GGA
    XC_ABGGA   =64, // Generic GGA alpha,beta variables 
    XC_RSGGA   =128, // Generic GGA R,S variables
    XC_ABMGGA  =256, // Alpha/beta meta gga
    XC_RSMGGA  =512,
  };

#define XC_ALL_LDA (XC_ALDA | XC_RLDA | XC_ABLDA | XC_RSLDA)
#define XC_ALL_GGA (XC_AGGA | XC_RGGA | XC_ABGGA | XC_RSGGA)
#define XC_ALL_MGGA (XC_ABMGGA | XC_RSMGGA)
#define XC_MODES_LEN 9

enum xc_functionals
  {
    XC_VWN5_CORRELATION = 0,
    XC_SLATER_EXCHANGE,
    XC_LYP_CORRELATION,
    XC_BECKE_EXCHANGE, //Not including the Slater part.
    XC_PBE_CORRELATION,
    XC_PBE_EXCHANGE,
    XC_PW92_CORRELATION,
    XC_SRLDA_ERF_EXCHANGE,
    XC_SRLDA_ERF_CORRELATION,
    /*	  XC_PW91_EXCHANGE, */
    /*
    
    XC_REVPBE_EXCHANGE_UNTESTED,
    KIN_TF_UNTESTED,       // (orbital-free) kinetic energy functionals next
    KIN_PW91_UNTESTED,*/
    XC_FUNLIST_LEN // This has to be last!
  };

int xc_parse_functional(const char *functional_name);

// Set the functional as defined by the parameters in the
// params array. 
void xc_set_functional(enum xc_mode mode, 
		       const double params[XC_FUNLIST_LEN]);

//Return the total number of terms calculated by xc_eval
int xc_len(enum xc_mode mode, int order); 

// Evaluate the exchange-correlation functional defined by
// a previous call to xc_set_functional() with the density
// parameters in densvar, with derivatives up to the given
// order. The result is put in the result array, in graded
// lexicographical ordering. This function also has a fortran 
// version.
void xc_eval(double result[], int order, const double densvars[]);

//Return the index (into result) for a given derivative, i.e.
//{0,2,1,0} gives the d^3/dS^2dZ term index.
int xc_index(enum xc_mode mode, const int derivative[]);

#ifdef WITH_QD
#include <qd/qd_real.h>
void xc_eval_qd(qd_real *result, int order, const qd_real densvars[]);
#endif

#endif
