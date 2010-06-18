#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xcfun.h"

#define TEST_ORDER 2

int main(int argc, char *argv[])
{
  xc_functional fun = xc_new_functional();
  if (argc > 1)
    {
      int i,j;
      for (i=1;i+1<argc;i+=2)
	{
	  double w;
	  int param = -1;
	  // Find a setting with name argv[i]
	  for (j=0;j<XC_NR_PARAMS;j++)
	    if (xc_name(j) && strcmp(xc_name(j),argv[i])==0)
	      {
		param = j;
		break;
	      }
	  if (param == -1)
	    {
	      fprintf(stderr,"Invalid setting '%s', quitting.\n",argv[i]);
	      return EXIT_FAILURE;
	    }
	  // Parse the value of the setting
	  if (sscanf(argv[i+1],"%lf",&w) != 1)
	    {
	      fprintf(stderr,
		      "Error parsing weight '%s', quitting.\n",argv[i+1]);
	      return EXIT_FAILURE;
	    }
	  // Set it
	  xc_set_param(fun,param,w);
	}
    }
  else
    {
      int i = 0;
      printf("Available functionals and settings:\n");
      for (i=0;i<XC_NR_PARAMS;i++)
	printf("%s\n",xc_name(i));
      printf("Usage: funeval FUNCTIONAL WEIGHT [FUNCTIONAL WEIGHT ..]\n");
      return 0;
    }
  xc_set_mode(fun,XC_VARS_AB);
  printf("Type is %i (0: LDA, 1: GGA, 2: MetaGGA)\n",xc_get_type(fun));
  printf("Maximal compiled order: %i\n", xc_max_order(fun));
  printf("Output length at order %i: %i\n",TEST_ORDER, 
	 xc_output_length(fun,TEST_ORDER));
  printf("Reading input density.. (%i values)\n",xc_input_length(fun));
  if (xc_max_order(fun) >= TEST_ORDER)
    {
      int i;
      double *inp = malloc(sizeof*inp*xc_input_length(fun));
      double *out = malloc(sizeof*out*xc_output_length(fun,TEST_ORDER));
      for (i=0;i<xc_input_length(fun);i++)
	if (scanf("%lf",&inp[i]) != 1)
	  {
	    fprintf(stderr,"Error reading density value, quitting.\n");
	    return EXIT_FAILURE;
	  }
      xc_eval(fun,TEST_ORDER,1,inp,out);
      for (i=0;i<xc_output_length(fun,TEST_ORDER);i++)
	printf("%f\n",out[i]);
    }
  return 0;
}
