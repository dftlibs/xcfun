#include <iostream>

using namespace std;

#include "taylor.hpp"
typedef double num_t;

using namespace polymul;

// Some taylor coefficients around x = 2.5, calculated using
// Maxima bfloat(taylor(...))
const num_t exp_good[] = {1.218249396070347e1,
                          1.218249396070347e1,
                          6.091246980351737e0,
                          2.030415660117246e0,
                          5.076039150293114e-1,
                          1.015207830058623e-1,
                          1.692013050097705e-2};

const num_t log_good[] = {9.162907318741551e-1,
                          4.0e-1,
                          -8.0e-2,
                          2.133333333333333e-2,
                          -6.4e-3,
                          2.048e-3,
                          -6.826666666666667e-4};

const num_t atan_good[] = {1.1902899496825317,
                           1.379310344827586e-1,
                           -4.756242568370987e-2,
                           1.552612516571679e-2,
                           -4.750587107528691e-3,
                           1.336092873197889e-3,
                           -3.310338712605161e-4};

// exp(1/(2.5 + 3x + 7y))
const num_t composed_good[] = {1.49182469764127,
                               -7.160758548678098e-1, // x
                               -1.670843661358223,    // y
                               1.031149231009646,     // x^2
                               4.812029744711682,     // xy
                               5.614034702163628};    // y^2

template <class num, int Nvar, int Ndeg>
void printpoly(ostream & dst, const polynomial<num, Nvar, Ndeg> & p) {
  int exps[Nvar] = {0};
  for (int i = 0; i < p.size; i++) {
    p.exponents(i, exps);
    for (int j = 0; j < Nvar; j++)
      cout << exps[j] << " ";
    cout << " " << p[i] << endl;
    assert(i == p.term_index(exps));
  }
}

#define NR_COEFF_CHECK 6

template <class T, int Nvar, int Ndeg>
int taylor_check(const char * label,
                 const taylor<T, Nvar, Ndeg> & t,
                 const num_t correct[],
                 num_t thres) {
  int nfail = 0;
  for (int i = 0; i < NR_COEFF_CHECK; i++) {
    T err = 2 * fabs(t[i] - correct[i]) / (1 + fabs(correct[i]));
    if (err > thres) {
      cout << "Error in coefficient " << i << " of " << label
           << ". Correct:" << correct[i] << ", taylor:" << t[i]
           << " error: " << correct[i] - t[i] << endl;
      nfail++;
    }
  }
  return nfail;
}

template <class T> T error_measure(const T & x1, const T & x2) {
  return 2 * fabs(x1 - x2) / (1 + 0.5 * fabs(x1 + x2));
}

template <class T, int Nvar, int Ndeg>
int taylor_compare(const taylor<T, Nvar, Ndeg> & t1,
                   const taylor<T, Nvar, Ndeg> & t2,
                   num_t thres) {
  int nfail = 0;
  for (int i = 0; i < t1.size; i++) {
    T err = 2 * fabs(t1[i] - t2[i]) / (1 + 0.5 * fabs(t1[i] + t2[i]));
    if (err > thres) {
      cout << "Difference in coefficient " << i << ". Correct:" << t1[i]
           << ", taylor:" << t2[i] << " difference: " << t1[i] - t2[i] << endl;
      nfail++;
    }
  }
  return nfail;
}

template <class T, int Nvar, int Ndeg>
T taylor_sumsq(const taylor<T, Nvar, Ndeg> & t) {
  T s = 0;
  for (int i = 0; i < t.size; i++)
    s += t[i] * t[i];
  return s;
}

long fac(int n) {
  long f = 1;
  while (n > 1)
    f *= n--;
  return f;
}

long stupid_binomial(int n, int k) {
  if (k < 0 || k > n)
    return 0;
  else
    return fac(n) / (fac(k) * fac(n - k));
}

template <class T, int N> taylor<T, 2, N> binomial_generating_function(void) {
  taylor<T, 2, N> x(0, 0), y(0, 1);
  return 1 / (1 - (1 + x) * y); // generating function for the binomial coefficients
}

int main(void) {
  num_t x0 = 2.5;
  taylor<num_t, 1, 6> tin(x0, 0);
  int res = 0;

  cout.precision(16);

  // Test some elementary functions (also test composition etc)
  res += taylor_check("exp", exp(tin), exp_good, 1e-15);
  res += taylor_check("log", log(tin), log_good, 1e-15);
  res += taylor_check("atan", atan(tin), atan_good, 1e-15);

  // Test multidimensional multiply consistency
  taylor<num_t, 4, 5> tmul(0), tacc(1), t1, t2, t3;
  for (int i = 0; i < t1.size; i++) {
    t1[i] = i + 1;
    t2[i] = -3 + i * i;
    t3[i] = i + i * i + 1;
  }
  // multiplication should commute
  tacc = t1;
  tacc *= t2;
  tacc *= t3;

  tmul = t2;
  tmul *= t1;
  tmul *= t3;
  res += taylor_compare(tacc, tmul, 1e-16);
  // add and subtract a bit
  if (taylor_sumsq(tmul - tacc) > 1e-16) {
    cout << "Subtraction error\n";
    res++;
  }
  // Evaluate multivariate polynomial and check the results
  {
    taylor<num_t, 2, 2> p;
    for (int i = 0; i < p.size; i++)
      p[i] = i + 1;
    num_t x[2] = {3, 7.1};
    if (fabs(p.eval(x) - 473.26) / 473.26 > 1e-15) {
      cout << "Error evaluating multivariate, error " << p.eval(x) - 473.26 << endl;
      res++;
    }
    taylor<num_t, 4, 5> p1, p2, pp;
    for (int i = 0; i < 10; i++)
      p1[i] = i + 1;
    for (int i = 0; i < 10; i++)
      p2[i] = i + 1 + (i & 3);
    pp = p1 * p2;
    num_t r[4] = {3.1, 4, 5, 6};
    if (fabs(pp.eval(r) - p1.eval(r) * p2.eval(r)) > 1e-15 * fabs(pp.eval(r))) {
      cout << "Error evaluating multivariate:" << pp.eval(r) << " "
           << p1.eval(r) * p2.eval(r) << endl;
      res++;
    }
    // Evaluate with taylor vars
    taylor<num_t, 1, 2> pv[4], pr;
    for (int i = 0; i < 4; i++)
      pv[i][1] = i;
    pr = tin.eval(pv);
  }

  // composition
  taylor<num_t, 2, 2> p(2.5);
  p[1] = 3;
  p[2] = 7; // Now p = 2.5 + 3*x + 7*y
  taylor<num_t, 2, 2> ret = exp(1 / p);
  res += taylor_check("exp(1/(2.5+3*x+7*y))", ret, composed_good, 1e-15);

  // Term indexing
  taylor<num_t, 3, 3> index_tester;
  int exponents1[] = {2, 0, 1};
  int exponents2[] = {0, 1, 1};
  /* 1 x y z x^2 xy xz y^2 yz z^2 x^3 x^2y x^2z xy^2 xyz xz^2 y^3 y^2z yz^2 z^3
     0 1 2 3 4   5  6  7   8  9   10  11   12   13   14  15   16  17   18  19*/
  if (index_tester.term_index(exponents1) != 12 ||
      index_tester.term_index(exponents2) != 8) {
    cout << "Error in taylor::index_of()\n";
    res++;
  }

  // Generating function binomial coeffs.
  taylor<num_t, 2, 10> nchoosek = binomial_generating_function<num_t, 10>();
  for (int i = 0; i <= 10; i++)
    for (int j = 0; i + j <= 10; j++) {
      int ij[2] = {i, j};
      if (fabs(nchoosek[nchoosek.term_index(ij)] - stupid_binomial(j, i)) > 1e-15) {
        cout << "Error evaluating binomial coefficients by generating function\n";
        res++;
      }
    }
  // Shifting
  taylor<num_t, 1, 6> sinx(1e-3, 0), sint = sin(sinx);
  taylor<num_t, 1, 12> shift_sinx(0, 0), shift_sint = sin(shift_sinx);
  taylor<num_t, 1, 6> shifted;
  num_t dx = 1e-3;
  shift_sint.shift(shifted, &dx);
  res += taylor_compare(sint, shifted, 1e-15);
  // Tensoring.
  // Test with a normal taylor series that completely contains
  // the tensor terms.
  tensored_taylor<3, double, 1, 2>::type tens_t = 0, tens_res;
  taylor<double, 3, 6> tens_cover = 0, tens_cover_out;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        tens_t[i][j][k] = 1 + i + 13 * j + 23 * k;
        int ijk[3] = {i, j, k};
        tens_cover[tens_cover.term_index(ijk)] = tens_t[i][j][k];
      }
  tens_res = 2 * exp(1 / log(tens_t + 1)) + atan(tens_t);
  tens_cover_out = 2 * exp(1 / log(tens_cover + 1)) + atan(tens_cover);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) {
        int ijk[3] = {i, j, k};
        if (error_measure(tens_res[i][j][k],
                          tens_cover_out[tens_cover_out.term_index(ijk)]) > 1e-14) {
          cout << "tensoring error: " << tens_res[i][j][k] << " vs "
               << tens_cover_out[tens_cover_out.term_index(ijk)] << endl;
          res++;
        }
      }
  // acos
  cout << scientific;
  cout.precision(16);
  taylor<num_t, 1, 6> acosx(1e-4, 0), acosout = acos(acosx);
  for (int i = 0; i < acosout.size; i++)
    cout << i << " " << acosout[i] << endl;

  // Test fast functions
  taylor<double, 2, 12> xin(2, 0), yin(1, 1), tslow, tfast;
  xin = xin + 3 * pow(xin, 2) + 2 * yin + pow(yin, 3);
  tslow = 1 / xin;
  taylor_reciprocal(tfast, xin);
  if (taylor_compare(tfast, tslow, 1e-16)) {
    cerr << "Fast reciprocal error." << endl;
    cout << "tslow: " << tslow << endl;
    cout << "tfast: " << tfast << endl;
    res++;
  }

  return res;
}
