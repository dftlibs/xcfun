#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "functional.h"

#ifndef M_PI // M_PI is not standard for some reason
#define M_PI 3.14159265358979323846
#endif

namespace xc_constants
{
  static const parameter CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);

  // PBE constants.
  static const parameter param_gamma = (1-log(2.0))/(M_PI*M_PI);
  static const parameter param_beta_pbe_paper = 0.066725;
  static const parameter param_beta_accurate  = 0.06672455060314922;
  static const parameter param_beta_gamma = param_beta_accurate/param_gamma;
}

#endif
