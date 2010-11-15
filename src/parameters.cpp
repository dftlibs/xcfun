#include "xcint.h"

static parameter_data data[] =
  {{"XC_RANGESEP_MU", "Range reparation parameter mu (1/r0)", 0.4}};


const parameter_data *xcint_param_lookup(enum xc_parameter p)
{
  return &data[p-XC_NR_FUNCTIONALS];
}
