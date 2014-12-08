#ifndef CONFIG_H
#define CONFIG_H

// Use #define XC_NO_REGULARIZATION to turn off
// checks and balances for physical densities.

// Enable functionals still in development (probably buggy)
#define XCFUN_IN_DEVELOPMENT

// Use #define XCFUN_REF_PW92C to use inaccurate constants in
// PW92C. This matches the reference implementation.

// Use inaccurate mu value in pbe exchange.
// #define XCFUN_ACCURATE_PBEX_MU

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
