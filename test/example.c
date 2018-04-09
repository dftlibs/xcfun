#include "xcfun.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void) {
  xc_functional fun = xc_new_functional();

  double d_elements[8] = {1, 2.1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};

  int nout, i;

  double * output;

  xc_set(fun, "blyp", 0.9);
  xc_set(fun, "pbec", 0.1);

  xc_eval_setup(fun, XC_A_B_AX_AY_AZ_BX_BY_BZ, XC_PARTIAL_DERIVATIVES, 1);

  nout = xc_output_length(fun);

  output = malloc(sizeof(*output) * nout);

  xc_eval(fun, d_elements, output);

  for (i = 0; i < nout; i++)
    printf("%.8e\n", output[i]);

  free(output);
  xc_free_functional(fun);
  return EXIT_SUCCESS;
}
