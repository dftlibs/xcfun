#include <iostream>
#include <cassert>
#include <cmath>
#define TAYLOR_CXXIO
#include "taylor.h"

/*
  Newton-Raphson style solvers.
  /ulfek
 */

using namespace std;

// For Ndeg == 1 this is Newton's method, Ndeg == 2 Halley's
template<class T, int Ndeg>
T householder(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &), const T &x0, const T &rhs, const T &xthres)
{
  assert(Ndeg>=1);
  taylor<T,1,Ndeg> x(x0,0), r;
  T d;
  for (int i=0;i<5;i++)
    {
      r = 1/(f(x)-rhs);
      d = r[Ndeg-1]/r[Ndeg];
      x[0] += d;
      cout << "x = " << x[0] << endl;
    }
  return x[0];
}

// Newton only, but avoids the 1/(f-rhs) so it can be used when f - rhs has no constant term
template<class T>
T newton(taylor<T,1,1> (*f)(const taylor<T,1,1> &), const T &x0, const T &rhs, const T &xthres)
{
  taylor<T,1,1> x(x0,0), r;
  T d;
  for (int i=0;i<5;i++)
    {
      r = f(x)-rhs;
      d = -r[0]/r[1];
      x[0] += d;
      cout << "x = " << x[0] << endl;
    }
  return x[0];
}


template<class T>
T f(const T &y)
{
  return exp(y);
}

int main()
{
  cout.precision(2);
  cout << scientific;
  taylor<double,1,17> rhs(1,0);
  newton<taylor<double,1,17> >(f,0.0,rhs,1e-14);
  return 0;
}
