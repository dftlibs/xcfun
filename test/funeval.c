#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xcfun.h"

#define TEST_ORDER 2

// Started with m = {0,0, .. 0} (nvar zeros) this function
// will increment the exponents in the m array to give the
// next derivative "exponents" for the coefficients in the 
// xc_eval output array.
static void next_exponents(int nvar, int m[])
{
  int i,k = 0;
  if (nvar == 1)
    {
      m[0]++;
      return;
    }
  for (i=0;i<nvar-1;i++)
    k += m[i];
  if (k == 0)
    {
      m[0] = m[nvar-1] + 1;
      m[nvar-1] = 0;
      return;
    }
  if (m[nvar-2] > 0)
    {
      m[nvar-1]++;
      m[nvar-2]--;
    }
  else
    {
      next_exponents(nvar-1,m);
      for (i=nvar-2;i>=0;i--)
	{
	  if (m[i] > 0)
	    {
	      m[i] += m[nvar-1];
	      break;
	    }
	}
      m[nvar-1] = 0;
    }
}

int main(int argc, char *argv[])
{
  xc_functional fun = xc_new_functional();
  int nvar, mode = XC_VARS_AB;
  int dobench = 0;
  printf("XCFun version: %g\n",xcfun_version());
  if (argc > 1)
    {
      int i,j;
      for (i=1;i+1<argc;i+=2)
	{
	  double w;
	  int param = -1;
	  if (strcmp(argv[i],"--benchmark") == 0)
	    {
	      dobench = 1;
	      i++;
	    }
	  if (strcmp(argv[i],"--ns") == 0)
	    {
	      mode = XC_VARS_NS;
	      i++;
	    }
	  if (strcmp(argv[i],"--n") == 0)
	    {
	      mode = XC_VARS_N;
	      i++;
	    }
	  if (strcmp(argv[i],"--ab") == 0)
	    {
	      mode = XC_VARS_AB;
	      i++;
	    }
	  // Find a setting with name argv[i]
	  for (j=0;j<XC_NR_PARAMS;j++)
	    if (xc_name(j) && (strcmp(xc_name(j),argv[i]) == 0))
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
  xc_set_mode(fun,mode);
  printf("Mode is %i (0: A, 1: N, 2: AB, 3: NS)\n",mode);
  printf("Type is %i (0: LDA, 1: GGA, 2: MetaGGA 3: MetaGGA[Tau+Laplacian])\n",xc_get_type(fun));
  printf("Maximal compiled order: %i\n", xc_max_order(fun));
  printf("Output length at order %i: %i\n",TEST_ORDER, 
	 xc_output_length(fun,TEST_ORDER));
  nvar = xc_input_length(fun);
  printf("Reading input density.. (%i values)\n",nvar);
  if (xc_max_order(fun) >= TEST_ORDER)
    {
      int i,j;
      double *inp = malloc(sizeof*inp*xc_input_length(fun));
      int *m = malloc(sizeof*m*xc_input_length(fun));
      double *out = malloc(sizeof*out*xc_output_length(fun,TEST_ORDER));
      for (i=0;i<xc_input_length(fun);i++)
	if (scanf("%lf",&inp[i]) != 1)
	  {
	    fprintf(stderr,"Error reading density value, quitting.\n");
	    return EXIT_FAILURE;
	  }
      // Only one point, so pitch is unimportant
      if (dobench)
	for (i = 1;i<1e6;i++)
	  xc_eval(fun,TEST_ORDER,inp,out); 
      else
	  xc_eval(fun,TEST_ORDER,inp,out); 
      for (i=0;i<nvar;i++)
	m[i] = 0;
      printf("Derivative        Value\n");
      for (i=0;i<xc_output_length(fun,TEST_ORDER);i++)
	{
	  for (j=0;j<nvar;j++)
	    printf("%i ",m[j]);
	  printf("  %.15f\n",out[i]);
	  next_exponents(nvar,m);
	}
    }
  else
    {
      fprintf(stderr,"ERROR: Derivative order too high, check config.h\n");
      fprintf(stderr,"max order is %i, you asked for %i\n",xc_max_order(fun),
	      TEST_ORDER);
      return EXIT_FAILURE;
    }
  return 0;
}
