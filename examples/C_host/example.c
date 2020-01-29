/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2020 Ulf Ekström and contributors.
 *
 * This file is part of XCFun.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * XCFun library, see: <https://xcfun.readthedocs.io/>
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "XCFun/xcfun.h"

int main(void) {
  xcfun_t * fun = xcfun_new();

  double d_elements[8] = {1, 2.1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};

  int nout, i;

  double * output;

  xcfun_set(fun, "blyp", 0.9);
  xcfun_set(fun, "pbec", 0.1);

  xcfun_eval_setup(fun, XC_A_B_AX_AY_AZ_BX_BY_BZ, XC_PARTIAL_DERIVATIVES, 1);

  nout = xcfun_output_length(fun);

  output = malloc(sizeof(*output) * nout);

  xcfun_eval(fun, d_elements, output);

  for (i = 0; i < nout; i++)
    printf("%.8e\n", output[i]);

  free(output);
  xcfun_delete(fun);
  return EXIT_SUCCESS;
}
