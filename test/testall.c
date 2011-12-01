#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "xcfun.h"

/*
  Run all tests for all functionals.
  This program is in C to test the no stdc++ feature of xcfun.
 */

void check(const char *what, int cond)
{
  if (!cond)
    {
      fprintf(stderr,"Failed check: %s\n",what);
      exit(-1);
    }
}

void checknum(const char *what, double x, double xref, double abserr, double relerr)
{
  if (fabs(x - xref) < abserr)
    return;
  if (fabs(x - xref) < relerr*xref)
    return; 
  fprintf(stderr,"Numerical check failed: %s\n got      %.16e\n expected %.16e\n Absolute error is %.4e, relative error is %.4e\n",
	  what,x,xref,x - xref,(x - xref)/xref);
  exit(-1);
}

/* Test permutation symmetries over variables and modes etc. */
void consistency_test(void)
{
  xc_functional fun = xc_new_functional();
  double d_unpolarized[8] = {1,1, 2,-3,4, 2,-3,4};
  double d_pol_a[8] = {1,2.1, 2,-3,4, 7,-8,9};
  double d_pol_b[8] = {2.1,1, 7,-8,9, 2,-3,4};
  int nout,i;
  double *output;
  double *out2;
  xc_set(fun,"pbe",1.0);
  xc_eval_setup(fun,XC_A_B_AX_AY_AZ_BX_BY_BZ,XC_PARTIAL_DERIVATIVES,1);
  nout = xc_output_length(fun);
  printf("%i\n",nout);
  check("correct output length 1",nout == 9); // 1 + 8
  output = malloc(sizeof(*output)*nout);
  out2 = malloc(sizeof(*output)*nout);
  xc_eval(fun,d_unpolarized,output);
  for (i=0;i<nout;i++)
    printf("%.14e\n",output[i]);
  checknum("unpolarized symmetry 1",output[1] - output[2],0,1e-14,1e-12);
  checknum("unpolarized symmetry 2",output[3] - output[6],0,1e-14,1e-12);
  checknum("unpolarized symmetry 3",output[4] - output[7],0,1e-14,1e-12);
  checknum("unpolarized symmetry 4",output[5] - output[8],0,1e-14,1e-12);
  xc_eval(fun,d_pol_a,output);
  xc_eval(fun,d_pol_b,out2);
  checknum("polarized symmetry 1",output[1] - out2[2],0,1e-14,1e-12);
  checknum("polarized symmetry 2",output[3] - out2[6],0,1e-14,1e-12);
  checknum("polarized symmetry 3",output[4] - out2[7],0,1e-14,1e-12);
  checknum("polarized symmetry 4",output[5] - out2[8],0,1e-14,1e-12);
  checknum("polarized symmetry 5",out2[1] - output[2],0,1e-14,1e-12);
  checknum("polarized symmetry 6",out2[3] - output[6],0,1e-14,1e-12);
  checknum("polarized symmetry 7",out2[4] - output[7],0,1e-14,1e-12);
  checknum("polarized symmetry 8",out2[5] - output[8],0,1e-14,1e-12);
  free(output);
  free(out2);

}

int main(void)
{
  int i = 0;
  const char *n,*s;
  consistency_test();
  printf("%s",xcfun_splash());
  printf("XCFun version: %g\n",xcfun_version());
  printf("\nAvailable functionals and other settings:\n");
  while ((n = xc_enumerate_parameters(i++)))
    {
      printf("%s \t",n);
      if ((s = xc_describe_short(n)))
	printf("%s",s);
      else
	printf("[No description]");
      printf("\n");
    }
  printf("\nAvailable aliases:\n");
  i = 0;
  while ((n = xc_enumerate_aliases(i++)))
    {
      printf("%s \t",n);
      if ((s = xc_describe_short(n)))
	printf("%s",s);
      else
	printf("[No description]");
      printf("\n");
    }
  printf("\nRunning tests..\n");
  if (xcfun_test() == 0)
    {
      printf("\nAll tests ok\n");
      return 0;
    }
  else
    {
      printf("\nSome tests failed\n");
      return -1;
    }
}
