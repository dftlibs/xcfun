#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xcfun.h"

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
  int dobench = 0;
  int quiet = 0;
  int order = 1;
  int test_res = xcfun_test();
  printf("nr failed tests: %i\n",test_res);
#if 0
  if (argc > 1)
    {
      int i,j;
      for (i=1;i<argc;i++)
	{
	  double w;
	  int param = -1;
	  if (strcmp(argv[i],"--benchmark") == 0)
	    {
	      dobench = 1;
	    }
	  else if (strcmp(argv[i],"--quiet") == 0)
	    {
	      quiet = 1;
	    }
	  else if (strcmp(argv[i],"--order") == 0)
	    {
	      if (!argv[i+1])
		{
		  fprintf(stderr,"Expected integer after --order, quitting.\n");
		  return EXIT_FAILURE;
		}
	      if (sscanf(argv[i+1],"%i",&order) != 1)
		{
		  fprintf(stderr,"Error reading the order from '%s', quitting.\n",argv[i+1]);
		  return EXIT_FAILURE;
		}
	      i++;
	    }
	  // Find a setting with name argv[i]
	  else
	    {
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
	      if (!argv[i+1])
		{
		  fprintf(stderr,"Expected number after '%s', quitting.\n",argv[i]);
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
	      i++;
	    }
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
#endif
  xc_set(fun,XC_LYPC,1.0);
  xc_set_mode(fun,XC_PARTIAL_DERIVATIVES);
  if (!xc_try_vars(fun,XC_A_B_GAA_GAB_GBB))
    {
      fprintf(stderr,"Cannot set vars, quitting\n");
      return EXIT_FAILURE;
    }
  if (!xc_try_order(fun,2))
    {
      fprintf(stderr,"Cannot set order, quitting\n");
      return EXIT_FAILURE;
    }
#ifdef OLD
  nvar = xc_input_length(fun);
  if (!quiet)
    {
      printf("XCFun version: %g\n",xcfun_version());
      printf("Mode is %i (0: A, 1: N, 2: AB, 3: NS)\n",mode);
      printf("Type is %i (0: LDA, 1: GGA, 2: MetaGGA 3: MetaGGA[Tau+Laplacian])\n",xc_get_type(fun));
      printf("Maximal compiled order: %i\n", xc_max_order(fun));
      printf("Output length at order %i: %i\n",order, 
	     xc_output_length(fun,order));
      printf("Reading input density.. (%i values)\n",nvar);
    }
#endif
  int i,j;
  int inlen = 5;
  int outlen = 21;
  double *inp = malloc(sizeof*inp*inlen);
  int *m = malloc(sizeof*m*outlen);
  double *out = malloc(sizeof*out*outlen);
  for (i=0;i<inlen;i++)
    if (scanf("%lf",&inp[i]) != 1)
      {
	fprintf(stderr,"Error reading density value, quitting.\n");
	return EXIT_FAILURE;
      }
  // Only one point, so pitch is unimportant
  if (dobench)
    for (i = 1;i<1e6;i++)
      xc_eval(fun,inp,out); 
  else
    xc_eval(fun,inp,out); 
  for (i=0;i<inlen;i++)
    m[i] = 0;
  if (!quiet) 
    {
      printf("Derivative        Value\n");
      for (i=0;i<outlen;i++)
	{
	  for (j=0;j<inlen;j++)
	    printf("%i ",m[j]);
	  printf("  %.15f\n",out[i]);
	  next_exponents(inlen,m);
	}
    }
  else
    {
      for (i=0;i<outlen;i++)
	printf("%.16e ",out[i]);
      printf("\n");
    }
  return EXIT_SUCCESS;
}
