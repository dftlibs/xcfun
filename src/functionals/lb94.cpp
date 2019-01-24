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

#if 0 // Does not work and should maybe not be here in the first place
#include "functional.hpp"
#include "slater.hpp"
#include "vwn.hpp"

//Have to be a bit special here since lb94 is not a functional but 
//a potential.

template<typename T, int Nvar, int Ndeg>
static taylor<T,Nvar,Ndeg> lb94(const densvars<taylor<T,Nvar,Ndeg> > &d)
{
  const parameter beta = 0.05;
  taylor<T,Nvar,Ndeg> res = d.n*vwn::vwn5_eps(d) + slaterx(d);
  // LB94 is basically LDA with a GGA modification _for the potential_
  if (Ndeg >= 1)
    {
      // Alpha potential
      T x = d.gaa[0]*pow(d.a[0],-4.0/3.0);
      T v = -beta*pow(d.a[0],1.0/3.0)*x*x/(1+3*beta*x*asinh(x)); 
      // Use the chain rule to calculate contributions to
      // each dE/dx for the unknown x. 
      for (int i=0;i<Nvar;i++)
	res[1+i] += v*d.a[1+i];

      // Beta potential
      x = d.gbb[0]*pow(d.b[0],-4.0/3.0);
      v = -beta*pow(d.b[0],1.0/3.0)*x*x/(1+3*beta*x*asinh(x)); 
      for (int i=0;i<Nvar;i++)
	res[1+i] += v*d.b[1+i];
    }
  return res;
}

void setup_lb94(functional &f)
{
  f.describe(XC_LB94,XC_GGA,
	     "LB94 Exchange-correlation functional",
	     "Exchange-correlation potential with correct asymtotic behavior\n"
	     "R. van Leeuwen and E. J. Baerends PRA 49, 2421 (1994)\n"
	     "Note that the LB94 energy is not well defined, here its"
	     "simply the SVWN5 energy, also used for order >= 2.\n"
	     "Implemented by Ulf Ekstrom.\n");

  SET_GGA_ENERGY_FUNCTION(f, lb94);

  const double d[5] =
    {
     0.39e+02,
     0.38e+02,
     0.81e+06,
     0.82e+06,
     0.82e+06
    };

  const double ref[21] =
    {
//     radovan: reference data obtained from *.c implementation in DIRAC
      -2.504589269450e+02, // 00000
      -4.326578426659e+00, // 10000
      -4.292112232065e+00, // 01000
       0.000000000000e+00, // 00100
       0.000000000000e+00, // 00010
       0.000000000000e+00, // 00001
      -3.520452674168e-02, // 20000
      -1.028611984896e-03, // 11000
       0.000000000000e+00, // 10100
       0.000000000000e+00, // 10010
       0.000000000000e+00, // 10001
      -3.578939264344e-02, // 02000
       0.000000000000e+00, // 01100
       0.000000000000e+00, // 01010
       0.000000000000e+00, // 01001
       0.000000000000e+00,
       0.000000000000e+00,
       0.000000000000e+00,
       0.000000000000e+00,
       0.000000000000e+00,
       0.000000000000e+00
    };

  f.add_test(XC_VARS_AB, 2, d, ref, 1e-11);
}
#endif
