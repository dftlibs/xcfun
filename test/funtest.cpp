#include <iostream>
#include <qd/qd_real.h>
#include "taylor.h"
#include "xc_fun.h"

using namespace std;

template<class num, int Nvar, int Ndeg>
void print_rms_diff(const polynomial<num,Nvar,Ndeg> &p1,
                    const polynomial<num,Nvar,Ndeg> &p2)
{
  cout << "Order \tRMS \t\tRMS_Error\tRelative_Error" << endl;
  for (int n=0;n<=Ndeg;n++)
    {
      num rms = 0.0, avg = 0.0;
      for (int i=polylen(Nvar,n-1);i<polylen(Nvar,n);i++)
        {
          rms += pow(p1[i] - p2[i],2);
          avg += p1[i]*p1[i] + p2[i]*p2[i];
        }
      avg = sqrt(avg/(2*polylen(Nvar-1,n)));
      rms = sqrt(rms/polylen(Nvar-1,n));
      cout << n << " \t" << avg << " \t" << rms <<  " \t" << rms/avg << endl;
    }
}

template<class num, int Nvar, int Ndeg>
void print_rms_diff_compact(const polynomial<num,Nvar,Ndeg> &p1,
                    const polynomial<num,Nvar,Ndeg> &p2)
{
  for (int n=0;n<=Ndeg;n++)
    {
      num rms = 0.0, avg = 0.0;
      for (int i=polylen(Nvar,n-1);i<polylen(Nvar,n);i++)
        {
          rms += pow(p1[i] - p2[i],2);
          avg += p1[i]*p1[i] + p2[i]*p2[i];
        }
      avg = sqrt(avg/(2*polylen(Nvar-1,n)));
      rms = sqrt(rms/polylen(Nvar-1,n));
      cout << rms/avg << " ";
    }
  cout << endl;
}


#define TEST_ORDER 4

int main(void)
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

  fun[XC_VWN5_CORRELATION] = 1;
  //  fun[XC_BECKE_EXCHANGE] = 1;
  xc_set_functional(XC_ABFGH,fun);
  while (1)
    {
      cin >> dd[0] >> dd[1] >> dd[2] >> dd[3] >> dd[4];
#ifdef WITH_QD
      for (int i=0;i<5;i++)
	dqd[i] = dd[i];
#endif
      xc_eval(res_d.c,TEST_ORDER,dd);
#ifdef WITH_QD
      xc_eval(res_qd.c,TEST_ORDER,dqd);
      cout.precision(20);
      res_d.convert_to(cmp);
      print_rms_diff_compact(cmp,res_qd);
      fpu_fix_end(&cw);
#endif
    }
  return 0;
}
