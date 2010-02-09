#include "xc_fun.h"
#include "taylor.h"

#define ORDER 2

int main(void)
{
  double params[XC_PARAMS_LEN] = {0};
  params[XC_LYP_CORRELATION] = 1;
  params[XC_BECKE_EXCHANGE] = 1;
  xc_set_functional(params);
  double d[5] = {1,0.5,0.2,0.1,0.1};
  double res[taylor<double,5,ORDER>::size] = {0};
  for (int i=0;i<5e5;i++)
    xc_eval_gga(res,ORDER,d);
  return res[0];
}
