#include <iostream>
#include "xcfun.h"

using namespace std;

int main(void)
{
  xc_functional fun;
  int n = 0;
  const char *s;
  cout << "Available functionals and settings:" << endl;
  while ((s = fun.setting_name(n)))
    {
      cout << s << endl;
      n++;
    }
  cout << "Running tests.." << endl;
  if (xcfun_test() == 0)
    {
      cout << "All tests ok" << endl;
      return 0;
    }
  else
    {
      cout << "Some test(s) failed" << endl;
      return -1;
    }
}



