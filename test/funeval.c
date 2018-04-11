#include "xcfun.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Started with m = {0,0, .. 0} (nvar zeros) this function
// will increment the exponents in the m array to give the
// next derivative "exponents" for the coefficients in the
// xc_eval output array.
static void next_exponents(int nvar, int m[]) {
  int i, k = 0;
  if (nvar == 1) {
    m[0]++;
    return;
  }
  for (i = 0; i < nvar - 1; i++)
    k += m[i];
  if (k == 0) {
    m[0] = m[nvar - 1] + 1;
    m[nvar - 1] = 0;
    return;
  }
  if (m[nvar - 2] > 0) {
    m[nvar - 1]++;
    m[nvar - 2]--;
  } else {
    next_exponents(nvar - 1, m);
    for (i = nvar - 2; i >= 0; i--) {
      if (m[i] > 0) {
        m[i] += m[nvar - 1];
        break;
      }
    }
    m[nvar - 1] = 0;
  }
}

int main(int argc, char * argv[]) {
  xc_functional fun = xc_new_functional();
  int nvar, mode = XC_A_B;
  int quiet = 0;
  int order = 1;
  int rep_times = 1;
  if (argc > 1) {
    int i;
    for (i = 1; i < argc; i++) {
      double w;
      if (strcmp(argv[i], "--quiet") == 0) {
        quiet = 1;
      } else if (strcmp(argv[i], "--ns") == 0) {
        mode = XC_N_S;
      } else if (strcmp(argv[i], "--n") == 0) {
        mode = XC_N;
      } else if (strcmp(argv[i], "--ngnn") == 0) {
        mode = XC_N_GNN;
      } else if (strcmp(argv[i], "--ngnntaun") == 0) {
        mode = XC_N_GNN_TAUN;
      } else if (strcmp(argv[i], "--ab") == 0) {
        mode = XC_A_B;
      } else if (strcmp(argv[i], "--order") == 0) {
        if (!argv[i + 1]) {
          fprintf(stderr, "Expected integer after --order, quitting.\n");
          return EXIT_FAILURE;
        }
        if (sscanf(argv[i + 1], "%i", &order) != 1) {
          fprintf(
              stderr, "Error reading the order from '%s', quitting.\n", argv[i + 1]);
          return EXIT_FAILURE;
        }
        i++;
      } else if (strcmp(argv[i], "--rep") == 0) {
        if (!argv[i + 1]) {
          fprintf(stderr, "Expected integer after --rep, quitting.\n");
          return EXIT_FAILURE;
        }
        if (sscanf(argv[i + 1], "%i", &rep_times) != 1) {
          fprintf(stderr,
                  "Error reading the number of repetitions from '%s', quitting.\n",
                  argv[i + 1]);
          return EXIT_FAILURE;
        }
        i++;
      }
      // Find a setting with name argv[i]
      else {
        if (!argv[i + 1]) {
          fprintf(stderr, "Expected number after '%s', quitting.\n", argv[i]);
          return EXIT_FAILURE;
        }
        // Parse the value of the setting
        if (sscanf(argv[i + 1], "%lf", &w) != 1) {
          fprintf(stderr, "Error parsing weight '%s', quitting.\n", argv[i + 1]);
          return EXIT_FAILURE;
        }
        // Set it
        xc_set(fun, argv[i], w);
        i++;
      }
    }
  } else {
    printf("Usage: funeval FUNCTIONAL WEIGHT [FUNCTIONAL WEIGHT ..]\n");
    return 0;
  }

  int res = xc_eval_setup(fun, mode, XC_PARTIAL_DERIVATIVES, order);
  if (res != 0) {
    fprintf(stderr, "Error in setup, code %i. Quitting.\n", res);
    return EXIT_FAILURE;
  }
  nvar = xc_input_length(fun);
  while (1) {
    if (!quiet) {
      printf("XCFun version: %g\n", xcfun_version());
      printf("Mode is %i\n", mode);
      printf("Output length at order %i: %i\n", order, xc_output_length(fun));
      printf("Reading input density.. (%i values)\n", nvar);
    }
    int i, j;
    double * inp = malloc(sizeof *inp * xc_input_length(fun));
    int * m = malloc(sizeof *m * xc_input_length(fun));
    double * out = malloc(sizeof *out * xc_output_length(fun));
    for (i = 0; i < xc_input_length(fun); i++)
      if (scanf("%lf", &inp[i]) != 1) {
        if (!quiet && !feof(stdin))
          fprintf(stderr, "Error reading density value, quitting.\n");
        return EXIT_FAILURE;
      }
    // Only one point, so pitch is unimportant
    for (i = 0; i < rep_times; i++)
      xc_eval(fun, inp, out);
    for (i = 0; i < nvar; i++)
      m[i] = 0;
    if (!quiet) {
      printf("Derivative        Value\n");
      for (i = 0; i < xc_output_length(fun); i++) {
        for (j = 0; j < nvar; j++)
          printf("%i ", m[j]);
        printf("  %.15f\n", out[i]);
        next_exponents(nvar, m);
      }
    } else {
      for (i = 0; i < xc_output_length(fun); i++)
        printf("%.16e ", out[i]);
      printf("\n");
    }
  }
  return 0;
}
