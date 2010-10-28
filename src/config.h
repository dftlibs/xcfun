#ifndef CONFIG_H
#define CONFIG_H

//Maximum derivative order. Lower orders can be generated
//for GGA's and MGGA's, to avoid huge code size.
#define XC_MAX_ORDER 2
#define XC_LDA_MAX_ORDER XC_MAX_ORDER
#define XC_GGA_MAX_ORDER XC_MAX_ORDER
#define XC_MGGA_MAX_ORDER XC_MAX_ORDER
#define XC_MLGGA_MAX_ORDER XC_MAX_ORDER
#define XC_CONTRACT_MAX_ORDER XC_MAX_ORDER

// Use #define XC_NO_REGULARIZATION to turn off
// checks and balances for physical densities.

// Unfortunately many tests are designed for gradient
// values that violate the Cauchy-Schwarz inequality,
// so turn off this check. It is anyway probably not needed.
#define XC_NO_SCHWARZ_REGULARIZATION

// Enable functionals still in development (probably buggy)
#define XCFUN_IN_DEVELOPMENT

// Use #define XCFUN_REF_PW92C to use inaccurate constants in
// PW92C. This matches the reference implementation.

//#define XCFUN_REF_PW92C

// Use #define XCFUN_VWN5_PBEC to use VWN5 as the LDA energy in PBEC.
// This is used in the ADF program (for historical reasons)
//#define XCFUN_VWN5_PBEC

// Some platforms do not have erf(), make erf() calls stop
// AT RUN TIME!
#define XCFUN_NO_ERF 

// This is the internal scalar type of the library, can be
// different from the external interface. 
#ifndef WITH_QD
typedef double ireal_t; 
#define INNER_TO_OUTER(INNER) INNER
#else
#include <qd/qd_real.h>
typedef qd_real ireal_t;
#define XCFUN_NUM_CONVERT // Must convert real types at i/o
#define INNER_TO_OUTER(INNER) to_double(INNER)
#endif

#endif
