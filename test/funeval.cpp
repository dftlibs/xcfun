#include <iostream>
#ifdef WITH_WD
#include <qd/qd_real.h>
#endif
#include <cstdlib>
#include <cstdio>
#include "taylor.h"
#include "xc_fun.h"

using namespace std;

template<class num, int Nvar, int Ndeg>
void printpoly(ostream &dst, const polynomial<num,Nvar,Ndeg> &p)
{
  int exps[Nvar] = {0};
  for (int i=0;i<p.size;i++)
    {
      p.exponents(i,exps);
      for (int j=0;j<Nvar;j++)
	cout << exps[j] << " ";
      cout << " " << p[i] << endl;
    }
}

#define TEST_ORDER 2

int main(int argc, char *argv[])
{
#ifdef WITH_QD
  unsigned int cw;
  fpu_fix_start(&cw);
  taylor<qd_real, 5, TEST_ORDER> res_qd, cmp;
  qd_real dqd[5] = {1.0,0.0,1.0,0.0,0.0};
#endif

  double dd[5] =   {1.0,0.0,1.0,0.0,0.0};
  taylor<double, 5, TEST_ORDER> res_d; 
  double fun[XC_PARAMS_LEN] = {0};

#ifdef WITH_SINGLE
  float fd[5] =   {1.0,0.0,1.0,0.0,0.0};
  taylor<float, 5, TEST_ORDER> res_f; 
#endif

  for (int i=1;i+1<argc;i+=2)
    {
      int p = xc_parse_functional(argv[i]);
      if (p<0)
	{
	  cerr << "Unknown functional " << argv[i] << ", quitting." << endl;
	  return EXIT_FAILURE;
	}
      double w;
      if (sscanf(argv[i+1],"%lf",&w) != 1)
	{
	  cerr << argv[i+1] << " is not a float, quitting." << endl;
	  return EXIT_FAILURE;
	}
      cerr << "Setting " << argv[i] << " to " << w << endl;
      fun[p] = w;
    }
  xc_set_functional(XC_ABFGH,fun);
  cerr << "Input density: " << endl;
  cout.precision(20);
  cin >> dd[0] >> dd[1] >> dd[2] >> dd[3] >> dd[4];
#ifdef WITH_QD
  for (int i=0;i<2;i++)
    dqd[i] = dd[i];
#endif
#ifdef WITH_SINGLE
  for (int i=0;i<2;i++)
    fd[i] = dd[i];
#endif
  xc_eval(res_d.c,TEST_ORDER,dd);
  cout << "XC (double):" << endl;
  printpoly(cout,res_d);
#ifdef WITH_SINGLE
  xc_eval_single(res_f.c,TEST_ORDER,fd);
  cout << "XC (single):" << endl;
  printpoly(cout,res_f);
#endif
#ifdef WITH_QD
  xc_eval_qd(res_qd.c,TEST_ORDER,dqd);
  cout << "XC (quad double):" << endl;
  printpoly(cout,res_qd);
  fpu_fix_end(&cw);
#endif
  return 0;
}
