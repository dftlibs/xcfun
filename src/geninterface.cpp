/* Program to generate a C header and Fortran module file
   with the current settings of XCFun. */

#include <cstdio>
#include <cstdlib>
#include "xcfun_internal.h"

void xcint_die(const char *message, int code)
{
  fprintf(stderr,"XCFun fatal error %i: ",code);
  fprintf(stderr,"%s",message);
  fprintf(stderr,"\n");
  exit(-1);
}


int main()
{
  /* C interface */
  FILE *of = fopen("include/xcfun_autogen.h","w");
  fprintf(of,"enum xcfun_parameters {\n");
  for (int i=0;i<XC_NR_PARAMS;i++)
    fprintf(of,"%s,\n",xc_name(i));
  fprintf(of,"XC_NR_PARAMS\n};\n");
  fclose(of);
  
  /* Fortran interface */
  of = fopen("fortran/xcfun_autogen.F90","w");
  fprintf(of,
	  "module xcfun\n"
	  "  implicit none\n");
  fprintf(of,"  integer, parameter :: XC_NR_PARAMS = %i\n",XC_NR_PARAMS);
  for (int i=0;i<XC_NR_PARAMS;i++)
    fprintf(of,"  integer, parameter :: %s = %i\n",
	    xc_name(i),i);
  fprintf(of,"end module\n");
  fclose(of);
  return 0;
}
