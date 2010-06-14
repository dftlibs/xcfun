/* A test to see if we can link and use xcfun from C without
   libstdc++ */
#include <stdio.h>
#include <stdlib.h>
#include "xcfun.h"

int main(void)
{
  printf("%s",xcfun_splash());
  return EXIT_SUCCESS;
}
