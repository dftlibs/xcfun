/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2019 Ulf Ekström and contributors.
 *
 * This file is part of XCFun.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * For information on the complete list of contributors to the
 * XCFun library, see: <https://xcfun.readthedocs.io/>
 */

/* Needs some work to work with ctaylor */
#include "functional.hpp"
#include "slater.hpp"
#include "taylor.hpp"

// This is the function we want to find roots for
template <typename T> T BR_z(const T & x) {
  return (x - 2) / x * exp(2.0 / 3.0 * x);
}

static double NR_step(double x, double z) {
  return (x * (3 * x * (exp(-2.0 / 3.0 * x) * z - 1) + 6)) / (x * (2 * x - 4) + 6);
}

// Return an x satisfying BR_z(x) = z
static double BR(double z) {
  double x0;
  if (z < -1e4)
    x0 = -2 / z;
  else if (z < -2)
    x0 = (sqrt(9 * z * z + 6 * z + 49) + 3 * z + 1) / 4;
  else if (z < 1)
    x0 = 2 * (z * exp(-4.0 / 3.0) + 1);
  else
    x0 = 3.0 / 2.0 * log(z) + 3.75 / (1.5 + log(z));
  for (int i = 0; i < 20; i++) {
    double xold = x0;
    x0 += NR_step(x0, z);
    if (fabs(xold - x0) < 1e-15 * (1 + x0))
      return x0;
  }
  fprintf(stderr, "BR: Not converged for z = %e\n", z);
  return x0;
}

// Obtain the Taylor expansion of x(y), which is the
// inverse of BR_y. Use linear method for simplicity.
template <typename T, int Ndeg>
void BR_taylor(const T & z0, taylor<T, 1, Ndeg> & t) {
  taylor<T, 1, Ndeg> f, d;
  t = 0;
  t[0] = BR(z0);
  t[1] = 1;
  f = BR_z(t);
  t[1] = 1 / f[1];
  // Linear method, for quadratic see i.e. Brent & Kung ~197x
  for (int i = 2; i <= Ndeg; i++) {
    f = BR_z(t);
    t[i] = -f[i] * t[1];
  }
}

/* This is a fully differentiable solver for Eq.(21) in
   Becke and Roussel, PRA 39, 1989. t is the _reciprocal_ right hand
   side value, x is returned.
 */
template <typename T, int Nvar>
static ctaylor<T, Nvar> BR(const ctaylor<T, Nvar> & t) {
  taylor<T, 1, Nvar> tmp;
  BR_taylor(t.c[0], tmp);

  ctaylor<T, Nvar> res = tmp[0];
  for (int i = 1; i <= Nvar; i++)
    res += tmp[i] * pow(t - t.c[0], i);
  return res;
}

template <typename num>
static num polarized(const num & na,
                     const num & gaa,
                     const num & lapa,
                     const num & taua,
                     const num & jpaa) // Becke tau here, no factor 1/2
{
  // The original BR article has a gamma constant, here this is put to 1.0
  num Q = (lapa - 2 * taua + (0.5 * gaa + 2 * jpaa) / na) / 6.0;
  num x = BR((1.0 / (2.0 / 3.0 * pow(M_PI, 2.0 / 3.0))) * Q * pow(na, -5.0 / 3.0));
  num b = cbrt(pow3(x) * exp(-x) / (8 * M_PI * na));
  return -(1 - (1 + 0.5 * x) * exp(-x)) / b; // FIXME: use expm1
}

template <typename num> static num brx(const densvars<num> & d) {
  return 0.5 * (d.a * polarized(d.a, d.gaa, d.lapa, 2 * d.taua, d.jpaa) +
                d.b * polarized(d.b, d.gbb, d.lapb, 2 * d.taub, d.jpbb));
}

template <typename num> static num brc(const densvars<num> & d) {
  parameter cab = 0.63, caa = 0.88;
  num UXa = polarized(d.a, d.gaa, d.lapa, 2 * d.taua, d.jpaa);
  num UXb = polarized(d.b, d.gbb, d.lapb, 2 * d.taub, d.jpbb);
  num zaa = abs(caa * (2.0 / UXa));
  num zbb = abs(caa * (2.0 / UXb));
  num zab = abs(cab * (1.0 / UXa + 1.0 / UXb));
  num ECopp = -0.8 * d.a * d.b * zab * zab * (1 - log(1 + zab) / zab);
  num ECaa = -0.01 * d.a * (2 * d.taua - (0.25 * d.gaa + d.jpaa) / d.a) *
             pow(zaa, 4) * (1 - 2 / zaa * log(1 + zaa / 2));
  num ECbb = -0.01 * d.b * (2 * d.taub - (0.25 * d.gbb + d.jpbb) / d.b) *
             pow(zbb, 4) * (1 - 2 / zbb * log(1 + zbb / 2));
  return ECopp + ECaa + ECbb;
}

template <typename num> static num brxc(const densvars<num> & d) {
  parameter cab = 0.63, caa = 0.88;
  num UXa = polarized(d.a, d.gaa, d.lapa, 2 * d.taua, d.jpaa);
  num UXb = polarized(d.b, d.gbb, d.lapb, 2 * d.taub, d.jpbb);
  num zaa = abs(caa * (2.0 / UXa));
  num zbb = abs(caa * (2.0 / UXb));
  num zab = abs(cab * (1.0 / UXa + 1.0 / UXb));
  num ECopp = -0.8 * d.a * d.b * zab * zab * (1 - log(1 + zab) / zab);
  num ECaa = -0.01 * d.a * (2 * d.taua - (0.25 * d.gaa + d.jpaa) / d.a) *
             pow(zaa, 4) * (1 - 2 / zaa * log(1 + zaa / 2));
  num ECbb = -0.01 * d.b * (2 * d.taub - (0.25 * d.gbb + d.jpbb) / d.b) *
             pow(zbb, 4) * (1 - 2 / zbb * log(1 + zbb / 2));
  return 0.5 * (UXa * d.a + UXb * d.b) + ECopp + ECaa + ECbb;
}

FUNCTIONAL(XC_BRX) = {
    "Becke-Roussells exchange with jp dependence",
    "See Becke, Canadian Journal of Chemistry, 1996, 74(6): 995-997"
    "Implemented by Ulf Ekstrom\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC | XC_LAPLACIAN | XC_JP,
    ENERGY_FUNCTION(brx)};

FUNCTIONAL(XC_BRC) = {
    "Becke-Roussells correlation with jp dependence",
    "See Becke, Canadian Journal of Chemistry, 1996, 74(6): 995-997"
    "Implemented by Ulf Ekstrom\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC | XC_LAPLACIAN | XC_JP,
    ENERGY_FUNCTION(brc)};

FUNCTIONAL(XC_BRXC) = {
    "Becke-Roussells correlation with jp dependence",
    "See Becke, Canadian Journal of Chemistry, 1996, 74(6): 995-997"
    "Implemented by Ulf Ekstrom\n",
    XC_DENSITY | XC_GRADIENT | XC_KINETIC | XC_LAPLACIAN | XC_JP,
    ENERGY_FUNCTION(brxc)};
