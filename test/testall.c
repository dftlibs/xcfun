#include <stdio.h>
#include <stdlib.h>
#include "xcfun.h"

/*
  Run all tests for all functionals.
  This program is in C to test the no stdc++ feature of xcfun.
 */

int main(void)
{
  int i = 0;
  const char *n,*s;
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
