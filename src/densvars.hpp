/*
 * XCFun, an arbitrary order exchange-correlation library
 * Copyright (C) 2020 Ulf Ekström and contributors.
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

#pragma once

#include "XCFunctional.hpp"
#include "config.hpp"

// When regularizing we shouldn't touch the higher order
// parts of the density, so we need this.
template <typename T, int N> void regularize(ctaylor<T, N> & x) {
  if (x < xcfun::XCFUN_TINY_DENSITY)
    x.set(0, xcfun::XCFUN_TINY_DENSITY);
}

template <typename T> static void regularize(T & x) {
  if (x < xcfun::XCFUN_TINY_DENSITY)
    x = xcfun::XCFUN_TINY_DENSITY;
}

// Variables for expressing functionals, these are redundant because
// different functionals have different needs.
// TODO: Make sure all variables are handled in the switch.
template <typename T> struct densvars {
  // Fills all density variables that can be filled from vars. Length of d
  // depends on vars.
  densvars(const XCFunctional * parent, const T * d) {
    this->parent = parent;
    switch (parent->vars) {
      case XC_A_GAA:
        gaa = d[1];
        gab = 0;
        gbb = 0;
        gnn = gaa;
        gss = gaa;
        gns = gaa;
      case XC_A:
        b = 0;
        n = a;
        regularize(n);
        s = n;
        break;
      case XC_A_B_GAA_GAB_GBB_TAUA_TAUB:
        taua = d[5];
        taub = d[6];
        tau = taua + taub;
      case XC_A_B_GAA_GAB_GBB:
        gaa = d[2];
        gab = d[3];
        gbb = d[4];
        gnn = gaa + 2 * gab + gbb;
        gss = gaa - 2 * gab + gbb;
        gns = gaa - gbb;
      case XC_A_B:
        a = d[0];
        regularize(a);
        b = d[1];
        regularize(b);
        n = a + b;
        s = a - b;
        break;
      case XC_N_S_GNN_GNS_GSS_TAUN_TAUS:
        taua = d[5] + d[6];
        taub = d[5] - d[6];
        tau = taua + taub;
      case XC_N_S_GNN_GNS_GSS:
        gnn = d[2];
        gns = d[3];
        gss = d[4];
        gaa = 0.25 * (gnn + 2 * gns + gss);
        gab = 0.25 * (gnn - gss);
        gbb = 0.25 * (gnn - 2 * gns + gss);
      case XC_N_S:
        n = d[0];
        regularize(n);
        s = d[1];
        a = n + s;
        regularize(a);
        b = n - s;
        regularize(b);
        break;
      case XC_N_GNN_TAUN:
        taua = d[2] / 2;
        taub = d[2] / 2;
        tau = d[2];
      case XC_N_GNN:
        gnn = d[1];
        gss = 0;
        gns = 0;
        gaa = 0.25 * gnn;
        gab = gaa;
        gbb = gaa;
      case XC_N:
        n = d[0];
        regularize(n);
        s = 0;
        a = 0.5 * n;
        b = a;
        break;
      case XC_N_2ND_TAYLOR:
        lapa = 0.5 * (d[4] + d[7] + d[9]);
        lapb = lapa;
      case XC_N_NX_NY_NZ_TAUN:
        taua = d[4] / 2;
        taub = d[4] / 2;
        tau = d[4];
      case XC_N_NX_NY_NZ:
        gnn = d[1] * d[1] + d[2] * d[2] + d[3] * d[3];
        gss = 0;
        gns = 0;
        gaa = 0.25 * gnn;
        gab = gaa;
        gbb = gaa;
        n = d[0];
        regularize(n);
        s = 0;
        a = 0.5 * n;
        b = a;
        break;
      case XC_A_B_2ND_TAYLOR: // a gax gay gaz haxx haxy haxz hayy hayz hazz
                              // 0 1   2   3   4    5    6    7    8    9
        lapa = d[4] + d[7] + d[9];
        lapb = d[14] + d[17] + d[19];
        a = d[0];
        regularize(a);
        b = d[10];
        regularize(b);
        gaa = d[1] * d[1] + d[2] * d[2] + d[3] * d[3];
        gab = d[1] * d[11] + d[2] * d[12] + d[3] * d[13];
        gbb = d[11] * d[11] + d[12] * d[12] + d[13] * d[13];
        gnn = gaa + 2 * gab + gbb;
        gss = gaa - 2 * gab + gbb;
        gns = gaa - gbb;
        n = a + b;
        s = a - b;
        break;
      case XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB:
        taua = d[8];
        taub = d[9];
        tau = d[8] + d[9];
      case XC_A_B_AX_AY_AZ_BX_BY_BZ:
        //    0 1  2  3  4  5  6  7
        a = d[0];
        regularize(a);
        b = d[1];
        regularize(b);
        gaa = d[2] * d[2] + d[3] * d[3] + d[4] * d[4];
        gab = d[2] * d[5] + d[3] * d[6] + d[4] * d[7];
        gbb = d[5] * d[5] + d[6] * d[6] + d[7] * d[7];
        gnn = gaa + 2 * gab + gbb;
        gss = gaa - 2 * gab + gbb;
        gns = gaa - gbb;
        n = a + b;
        s = a - b;
        break;
      case XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS:
        taua = 0.5 * (d[8] + d[9]);
        taub = 0.5 * (d[8] - d[9]);
        tau = d[8];
      case XC_N_S_NX_NY_NZ_SX_SY_SZ:
        //    0 1  2  3  4  5  6  7
        n = d[0];
        regularize(n);
        s = d[1];
        a = 0.5 * (n + s);
        regularize(a);
        b = 0.5 * (n - s);
        regularize(b);
        gnn = d[2] * d[2] + d[3] * d[3] + d[4] * d[4];
        gss = d[5] * d[5] + d[6] * d[6] + d[7] * d[7];
        gns = d[2] * d[5] + d[3] * d[6] + d[4] * d[7];
        gaa = 0.25 * (gnn + 2 * gns + gss);
        gab = 0.25 * (gnn - gss);
        gbb = 0.25 * (gnn - 2 * gns + gss);
        break;
      case XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB:
        jpaa = d[9];
        jpbb = d[10];
      case XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB:
        lapa = d[5];
        lapb = d[6];
        taua = d[7];
        taub = d[8];
        tau = taua + taub;
        gaa = d[2];
        gab = d[3];
        gbb = d[4];
        gnn = gaa + 2 * gab + gbb;
        gss = gaa - 2 * gab + gbb;
        gns = gaa - gbb;
        a = d[0];
        regularize(a);
        b = d[1];
        regularize(b);
        n = a + b;
        s = a - b;
        break;
      default:
        xcfun::die("Illegal/Not yet implemented vars value in densvars()",
                   parent->vars);
    }
    zeta = s / n;
    r_s = pow(3.0 / (n * 4.0 * M_PI), 1.0 / 3.0); // (3/4pi)^1/3*n^(-1/3) !check
    n_m13 = pow(n, -1.0 / 3.0);
    a_43 = pow(a, 4.0 / 3.0);
    b_43 = pow(b, 4.0 / 3.0);
  }

  const XCFunctional * parent{nullptr};
  double get_param(enum xc_parameter p) const { return parent->settings[p]; }

  T a{static_cast<T>(0)};
  T b{static_cast<T>(0)};
  T gaa{static_cast<T>(0)};
  T gab{static_cast<T>(0)};
  T gbb{static_cast<T>(0)};
  T n{static_cast<T>(0)};     /// na+nb
  T s{static_cast<T>(0)};     /// na - nb
  T gnn{static_cast<T>(0)};   /// (grad n) ^ 2
  T gns{static_cast<T>(0)};   /// (grad n).(grad s)
  T gss{static_cast<T>(0)};   /// (grad s) ^ 2
  T tau{static_cast<T>(0)};   /// Kinetic energy density.
  T taua{static_cast<T>(0)};  /// Alpha kinetic energy density.
  T taub{static_cast<T>(0)};  /// Beta kinetic energy density.
  T lapa{static_cast<T>(0)};  /// Alpha Laplacian density.
  T lapb{static_cast<T>(0)};  /// Beta Laplacian density.
  T zeta{static_cast<T>(0)};  /// s/n
  T r_s{static_cast<T>(0)};   /// (3/4pi)^1/3*n^(-1/3)
  T n_m13{static_cast<T>(0)}; /// pow(n,-1.0/3.0)
  T a_43{static_cast<T>(0)};
  T b_43{static_cast<T>(0)}; /// pow(a,4.0/3.0), pow(b,4.0/3.0)
  T jpaa{static_cast<T>(0)}; /// square of the alpha paramagnetic current vector.
  T jpbb{static_cast<T>(0)}; /// square of the beta paramagnetic current vector.
};
