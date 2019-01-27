#include "XCFun/xcfun.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
  Run all tests for all functionals.
  This program is in C to test the no stdc++ feature of xcfun.
 */

void check(const char * what, int cond) {
  if (!cond) {
    fprintf(stderr, "Failed check: %s\n", what);
    exit(-1);
  }
}

void checknum(const char * what,
              double x,
              double xref,
              double abserr,
              double relerr) {
  if (fabs(x - xref) < abserr)
    return;
  if (fabs(x - xref) < relerr * xref)
    return;
  fprintf(stderr,
          "Numerical check failed: %s\n got      %.16e\n expected %.16e\n Absolute "
          "error is %.4e, relative error is %.4e\n",
          what,
          x,
          xref,
          x - xref,
          (x - xref) / xref);
  exit(-1);
}

/* Test permutation symmetries over variables and modes etc. */
void consistency_test(void) {
  xc_functional fun = xc_new_functional();
  double d_unpolarized[8] = {1, 1, 2, -3, 4, 2, -3, 4};
  double d_pol_a[8] = {1, 2.1, 2, -3, 4, 7, -8, 9};
  double d_pol_b[8] = {2.1, 1, 7, -8, 9, 2, -3, 4};
  int nout;
  double * output;
  double * out2;
  xc_set(fun, "pbe", 1.0);
  xc_eval_setup(fun, XC_A_B_AX_AY_AZ_BX_BY_BZ, XC_PARTIAL_DERIVATIVES, 1);
  nout = xc_output_length(fun);
  check("correct output length 1", nout == 9); // 1 + 8
  output = malloc(sizeof(*output) * nout);
  out2 = malloc(sizeof(*output) * nout);
  xc_eval(fun, d_unpolarized, output);

  checknum("unpolarized symmetry 1", output[1] - output[2], 0, 1e-14, 1e-12);
  checknum("unpolarized symmetry 2", output[3] - output[6], 0, 1e-14, 1e-12);
  checknum("unpolarized symmetry 3", output[4] - output[7], 0, 1e-14, 1e-12);
  checknum("unpolarized symmetry 4", output[5] - output[8], 0, 1e-14, 1e-12);
  xc_eval(fun, d_pol_a, output);
  xc_eval(fun, d_pol_b, out2);
  checknum("polarized symmetry 1", output[1] - out2[2], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 2", output[3] - out2[6], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 3", output[4] - out2[7], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 4", output[5] - out2[8], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 5", out2[1] - output[2], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 6", out2[3] - output[6], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 7", out2[4] - output[7], 0, 1e-14, 1e-12);
  checknum("polarized symmetry 8", out2[5] - output[8], 0, 1e-14, 1e-12);
  free(output);
  free(out2);
  xc_free_functional(fun);
}

// Test that gradient square norma and gradient elements modes are consistent
void gradient_forms_test(void) {
  xc_functional fun = xc_new_functional();
  double d_elements[8] = {1, 2.1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
  double d_sqnorm[5] = {d_elements[0], d_elements[1]};
  int nout, i;
  double * output;
  double * out2;
  xc_set(fun, "blyp", 1.0);

  xc_eval_setup(fun, XC_A_B_AX_AY_AZ_BX_BY_BZ, XC_PARTIAL_DERIVATIVES, 1);
  nout = xc_output_length(fun);

  check("correct output length 1", nout == 9); // 1 + 8
  output = malloc(sizeof(*output) * nout);
  xc_eval(fun, d_elements, output);

  xc_eval_setup(fun, XC_A_B_GAA_GAB_GBB, XC_PARTIAL_DERIVATIVES, 1);
  nout = xc_output_length(fun);

  check("correct output length 1", nout == 6); // 1 + 5
  out2 = malloc(sizeof(*out2) * nout);
  d_sqnorm[2] = d_elements[2] * d_elements[2] + d_elements[3] * d_elements[3] +
                d_elements[4] * d_elements[4];
  d_sqnorm[3] = d_elements[2] * d_elements[5] + d_elements[3] * d_elements[6] +
                d_elements[4] * d_elements[7];
  d_sqnorm[4] = d_elements[5] * d_elements[5] + d_elements[6] * d_elements[6] +
                d_elements[7] * d_elements[7];
  xc_eval(fun, d_sqnorm, out2);

  checknum("Grad modes energy", output[0] - out2[0], 0, 1e-14, 1e-12);
  checknum("Grad modes density derivs alpha", output[1] - out2[1], 0, 1e-14, 1e-12);
  checknum("Grad modes density derivs beta", output[2] - out2[2], 0, 1e-14, 1e-12);
  // d/dg_ax = d/dgaa * g_ax + d/dgab * g_bx
  for (i = 0; i < 3; i++) {
    checknum("Grad modes density grad alpha",
             output[3 + i] -
                 (2 * out2[3] * d_elements[2 + i] + out2[4] * d_elements[5 + i]),
             0,
             1e-14,
             1e-12);
    checknum("Grad modes density grad beta",
             output[6 + i] -
                 (2 * out2[5] * d_elements[5 + i] + out2[4] * d_elements[2 + i]),
             0,
             1e-14,
             1e-12);
  }

  free(output);
  free(out2);
  xc_free_functional(fun);
}

void user_setup_test() {
    xc_functional fun1 = xc_new_functional();
    xc_functional fun2 = xc_new_functional();
    xc_functional fun3 = xc_new_functional();
    xc_set(fun1, "lda",  1.0);
    xc_set(fun2, "pbe",  1.0);
    xc_set(fun3, "m06l", 1.0);
    int rval1 = xc_user_eval_setup(fun1, 0, 0, 0, 1, 0, 0, 0, 0);
    int rval2 = xc_user_eval_setup(fun2, 1, 1, 1, 1, 0, 0, 0, 1);
    int rval3 = xc_user_eval_setup(fun3, 2, 2, 2, 1, 0, 1, 0, 1);
    check("Functional 1 correctly set up", rval1 == 0);
    check("Functional 2 correctly set up", rval2 == 0);
    check("Functional 3 correctly set up", rval3 == 0);
    xc_free_functional(fun1);
    xc_free_functional(fun2);
    xc_free_functional(fun3);
}

void xc_get_test() {
    xc_functional fun = xc_new_functional();
    xc_set(fun, "B3LYP", 1.0);

    double s, b, lyp, vwn, exx, kt, foo;
    check("SLATERX is a valid functional"   , xc_get(fun, "SLATERX", &s) == 0);
    check("BECKECORRX is a valid functional", xc_get(fun, "BECKECORRX", &b) == 0);
    check("LYPC is a valid functional"      , xc_get(fun, "LYPC", &lyp) == 0);
    check("VWN5C is a valid functional"     , xc_get(fun, "VWN5C", &vwn) == 0);
    check("EXX is a valid functional"       , xc_get(fun, "EXX", &exx) == 0);
    check("KTX is a valid functional"       , xc_get(fun, "KTX", &kt) == 0);
    check("FOO is NOT a valid functional"   , xc_get(fun, "FOO", &foo) != 0);

    checknum("B3LYP contains 80% SLATERX"   , s  , 0.80, 1.0e-14, 1.0e-12);
    checknum("B3LYP contains 72% BECKECORRX", b  , 0.72, 1.0e-14, 1.0e-12);
    checknum("B3LYP contains 81% LYPC"      , lyp, 0.81, 1.0e-14, 1.0e-12);
    checknum("B3LYP contains 19% VWN5C"     , vwn, 0.19, 1.0e-14, 1.0e-12);
    checknum("B3LYP contains 20% EXX"       , exx, 0.20, 1.0e-14, 1.0e-12);
    checknum("B3LYP contains  0% KTX"       , kt , 0.00, 1.0e-14, 1.0);

    xc_free_functional(fun);
}

int main(void) {
  int i = 0;
  const char *n, *s;
  consistency_test();
  gradient_forms_test();
  user_setup_test();
  xc_get_test();
  printf("%s", xcfun_splash());
  printf("XCFun version: %g\n", xcfun_version());
  printf("\nAvailable functionals and other settings:\n");
  while ((n = xc_enumerate_parameters(i++))) {
    printf("%s \t", n);
    if ((s = xc_describe_short(n)))
      printf("%s", s);
    else
      printf("[No description]");
    printf("\n");
  }
  printf("\nAvailable aliases:\n");
  i = 0;
  while ((n = xc_enumerate_aliases(i++))) {
    printf("%s \t", n);
    if ((s = xc_describe_short(n)))
      printf("%s", s);
    else
      printf("[No description]");
    printf("\n");
  }
  printf("\nRunning tests..\n");
  if (xcfun_test() == 0) {
    printf("\nAll tests ok\n");
    return 0;
  } else {
    printf("\nSome tests failed\n");
    return -1;
  }
}
