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
  int i = 0;
  printf("Available functionals and settings:\n");
  for (i=0;i<XC_NR_PARAMS;i++)
    {
      printf("%s\n",xc_name(i));
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
