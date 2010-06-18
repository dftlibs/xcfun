/* Program to generate a C header and Fortran module file
   with the current settings of XCFun. */

#include <cstdio>
#include <cstdlib>

enum xc_parameters
  {
#define PARAM(name) name
#include "list_of_parameters.h"
#undef PARAM
  };
static const char *param_symbols[XC_NR_PARAMS+1] =
  {
#define PARAM(name) #name
#include "list_of_parameters.h"
#undef PARAM
  };

int main()
{
  /* C interface */
  FILE *of = fopen("include/xcfun_autogen.h","w");
  fprintf(of,"enum xcfun_parameters {\n");
  for (int i=0;i<XC_NR_PARAMS;i++)
    fprintf(of,"%s,\n",param_symbols[i]);
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
	    param_symbols[i],i);
  fprintf(of,"end module\n");
  fclose(of);
  return 0;
}
