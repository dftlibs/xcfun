#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "xcfun.h"
#include "polymul.h"

using namespace std;

#define TEST_ORDER 2

int main(int argc, char *argv[])
{
  xc_functional fun;
  if (argc > 1)
    {
      for (int i=1;i+1<argc;i+=2)
	{
	  double w;
	  if (sscanf(argv[i+1],"%lf",&w) != 1)
	    {
	      cerr << "Error parsing weight '" << argv[i+1] << "'" << endl;
	      return -1;
	    }
	  if (fun.set_setting(argv[i],w) != 0)
	    {
	      cerr << "No such parameter '" << argv[i] << "'" << endl;
	      return -1;
	    }
	}
    }
  else
    {
      int n = 0;
      const char *s;
      cout << "Available functionals and settings:" << endl;
      while ((s = fun.setting_name(n)))
	{
	  cout << s << endl;
	  n++;
	}
      cout << "Usage: funeval FUNCTIONAL WEIGHT [FUNCTIONAL WEIGHT ..]" << endl;
      return 0;
    }
  fun.set_mode(XC_VARS_AB);
  cout << "Type is " << fun.get_type() << endl;
  cout << "Output length " << fun.output_length(TEST_ORDER) << endl;
  cout << "Maximal supported order " << fun.get_max_order() << endl;
  cout << "Input length " << fun.input_length() << endl;
  cout << "Reading input density.." << endl;
  cout.precision(15);
  if (fun.get_max_order() >= TEST_ORDER)
    {
      double dv[fun.input_length()];
      for (int i=0;i<fun.input_length();i++)
	cin >> dv[i];
      double res[fun.output_length(TEST_ORDER)];
      fun.eval(res,TEST_ORDER,dv);
      int expn[fun.input_length()];
      for (int i=0;i<fun.output_length(TEST_ORDER);i++)
	{
	  
	  cout << res[i] << endl;
	}
    }
  return 0;
}
