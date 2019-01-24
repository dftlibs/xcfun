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

#include "constants.hpp"
#include "functional.hpp"
#include "pw92eps.hpp"

// This is [(1+zeta)^(2/3) + (1-zeta)^(2/3)]/2, reorganized.
template <typename num> static num phi(const densvars<num> & d) {
  return pow(2.0, -1.0 / 3.0) * d.n_m13 * d.n_m13 * (sqrt(d.a_43) + sqrt(d.b_43));
}

template <typename num> static num energy(const densvars<num> & d) {
  using xc_constants::param_gamma;
  const parameter beta = 0.046;
  const parameter alpha = 1.8;
  //  const parameter omega = 4.5; // not needed if you use the fit
  num bg = beta / param_gamma;
  num eps = pw92eps::pw92eps(d);
  num u = phi(d);
  num u3 = pow3(u);
  // d2 is t^2
  num d2 = pow(1.0 / 12 * pow(3, 5.0 / 6.0) / pow(M_PI, -1.0 / 6), 2.0) * d.gnn /
           (u * u * pow(d.n, 7.0 / 3.0));
  num tt = pow(d2, 0.5); // this is t
  num v = tt * u * pow(d.r_s / 3.0, -1.0 / 6.0);
  num v3 = pow3(v);
  //
  // The term containing abs(zeta) gives problems with the
  // automatic derivatives
  // I tried to implement it in different ways. At the end I used
  // I fit to abs(z)^omega which is not containing the
  // absolute value of zeta (third version). I cheched this fit
  // against my own code implementing abs(z)^omega and
  // analytical derivatives of it. It seems that the fit works fine...
  // Nevertheless, if it is possible to fix the problems with the
  // automatic derivation of abs(z) it would be better to revert to
  // first version below.
  //
  //  FIRST VERSION
  // this is the only way i found to make abs(d.zeta)
  // num zz = d.zeta;
  // const parameter aa = 0.0;
  //  if (zz < aa)
  //    {
  //      num zz = -zz;
  //    }
  // num zw = exp(log(zz)*omega);
  //
  // SECOND VERSION
  //  num zz = d.zeta;
  //  const parameter aa = 0.0;
  //  num zw = 0.0;
  //  if (zz < aa)
  //    {
  //      num zw = exp(log(-zz)*omega); //pow(zz,omega);
  //    }
  //  else
  //    {
  //      num zw = exp(log(zz)*omega);
  //    }
  //
  // THIRD VERSION
  // I fitted abs(zeta)^omega with omega=4.5 with sum_i a_i*zeta^(2i)
  // note that with the values that I use this expression is always >0 in the
  // interval [-1,1]. A simpler fits are also
  // num zw = (0.679803 + 0.326987*d.zeta*d.zeta)*d.zeta*d.zeta*d.zeta*d.zeta;
  // num zw = (0.572797 + 0.635809*d.zeta*d.zeta -
  // 0.210351*d.zeta*d.zeta*d.zeta*d.zeta)*d.zeta*d.zeta*d.zeta*d.zeta;
  num zw = (0.462757 + 1.30129 * d.zeta * d.zeta -
            1.59546 * d.zeta * d.zeta * d.zeta * d.zeta +
            1.19635 * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta -
            0.36519 * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta * d.zeta *
                d.zeta) *
           d.zeta * d.zeta * d.zeta * d.zeta;
  //
  num ff = exp(-alpha * v3 * zw);
  num A = bg / expm1(-eps / (param_gamma * u3));
  num d2A = d2 * A;
  num H = param_gamma * u3 * log(1 + bg * d2 * (1 + d2A) / (1 + d2A * (1 + d2A)));
  return d.n * (eps + ff * H);
}

FUNCTIONAL(XC_ZVPBESOLC) = {"zvPBEsol correlation Functional",
                            "zvPBEsol correlation Functional\n"
                            "L.A. Constantin, E. Fabiano, F. Della Sala ,\n"
                            "      J. Chem. Phys. 137, 194105 (2012) .\n"
                            "Implemented by Eduardo Fabiano\n",
                            XC_DENSITY | XC_GRADIENT,
                            ENERGY_FUNCTION(energy)};
