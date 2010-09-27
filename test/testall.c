#include <stdio.h>
#include <stdlib.h>
#include "xcfun.h"

/*
  Run all tests for all functionals.
  This program is in C to test the no stdc++ feature of xcfun.
 */

int main(void)
{
  printf("%s",xcfun_splash());
  printf("XCFun version: %g\n",xcfun_version());
  int i = 0;
  printf("\nAvailable functionals and settings:\n");
  for (i=0;i<XC_NR_PARAMS;i++)
    {
      const char *s;
      printf("%s \t",xc_name(i));
      if ((s = xc_short_description(i)))
	printf("%s",s);
      else
	printf("[No description]");
      if (!xc_is_functional(i))
	printf(" [parameter or disabled functional]");
      printf("\n");
    }
  printf("Running tests..\n");
  if (xcfun_test() == 0)
    {
      printf("All tests ok\n");
      return 0;
    }
  else
    {
      printf("Some tests failed\n");
      return -1;
    }
}
