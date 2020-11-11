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

#include "config.hpp"
#include "constants.hpp"
#include "functional.hpp"
#include "pw92eps.hpp"

/*! Common code for SCAN and SCAN-like functionals
 *
 *  Implemented by James Furness,
 *
 *  SCAN    - J. Sun, A. Ruzsinszky, and J. P. Perdew, Phys. Rev. Lett. 115, 036402
 * (2015) rSCAN   - A. P. Bartok and J. R. Yates, J. Chem. Phys. 150, 161101 (2019)
 *  r2SCAN  - J Furness, A Kaplan, J Ning, J Perdew, & J Sun; J. Chem. Phys. Lett.;
 * Accepted (DOI: 10.1021/acs.jpclett.0c02405)
 *
 *  The r++SCAN (aka "rppSCAN") and r4SCAN functionals are relatives of r2SCAN. They
 * are not recommended for general use.
 *
 *  r++SCAN, r4SCAN - J Furness, A Kaplan, J Ning, J Perdew, & J Sun; in preparation.
 */

namespace SCAN_eps {
template <class num>
static num get_SCAN_Fx(const num,
                       const num,
                       const num,
                       const int,
                       const int,
                       const int);
template <class num>
static num SCAN_X_Fx(const num, const num, const parameter, const int, const int);
template <class num>
static num r2SCAN_C(const densvars<num> &, const int, const int, const int);
template <class num>
static num scan_ec0(const num,
                    const num,
                    const num,
                    const parameter,
                    const parameter,
                    const parameter);
template <class num>
static num lda_0(const num, const parameter, const parameter, const parameter);
template <class num>
static num scan_ec1(const num,
                    const num,
                    const num,
                    const parameter[8],
                    const parameter,
                    const parameter,
                    const parameter,
                    const parameter,
                    const int);
template <class num>
static void get_lsda1(const num, const num, const num, num &, num &);
template <class num>
static void gcor2(const parameter[6], const num, const num, num &, num &);

template <class num> static num fx_unif(const num & d) {
  return (-0.75 * pow(3 / PI, 1.0 / 3.0)) * pow(d, 4.0 / 3.0);
}

template <class num>
static num get_SCAN_Fx(const num d_n,
                       const num d_g,
                       const num d_tau,
                       const int IALPHA,
                       const int IINTERP,
                       const int IDELFX) {
  const parameter ETA = 1.0e-3;
  const parameter TAU_R = 1.0e-4;
  const parameter A_REG = 1.0e-3;

  num tauw = d_g / (8.0 * d_n);

  num tauUnif = 0.0;
  if (IALPHA == 1) {
    tauUnif = ((3.0 / 10.0) * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 5.0 / 3.0)) + TAU_R;
  } else {
    tauUnif = (0.3 * pow(3 * PI2, 2.0 / 3.0) * pow(d_n, 5.0 / 3.0));
  }

  num alpha = 0.0;
  if (IALPHA == 0) {
    // alpha (SCAN)
    if (abs(d_tau - tauw) > 1.0e-14) {
      alpha = (d_tau - tauw) / tauUnif;
    }

  } else if (IALPHA == 1) {
    // alpha' (rSCAN)
    alpha = (d_tau - tauw) / tauUnif;
    alpha = pow(alpha, 3) / (pow2(alpha) + A_REG);

  } else if (IALPHA == 2) {
    // \bar{alpha} (r2SCAN, r4SCAN)
    if (abs(d_tau - tauw) > 1.0e-14) {
      alpha = (d_tau - tauw) / (tauUnif + ETA * tauw);
    }

  } else {
    // Unknown IALPHA choice
    printf("ERROR: Unknown IALPHA %d\n", IALPHA);
  }

  num p = 0.0;
  if (abs(d_g) > 1.0e-16) {
    p = d_g /
        (4.0 * pow(3.0 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
  } else {
    p = 1e-16/
        (4.0 * pow(3.0 * PI2, 2.0 / 3.0) * pow(d_n, 8.0 / 3.0));
  }

  num Fx = SCAN_X_Fx(p, alpha, ETA, IINTERP, IDELFX);

  return Fx;
}

template <class num>
static num SCAN_X_Fx(const num p,
                     const num alpha,
                     const parameter ETA,
                     const int IINTERP,
                     const int IDELFX) {
  const parameter A1 = 4.9479;
  const parameter K1 = 0.065;
  const parameter K0 = 0.174;
  const parameter MU = 10.0 / 81.0;

  const parameter IE_PARAMS[8] = {1.0,
                                  -0.667,
                                  -0.4445555,
                                  -0.663086601049,
                                  1.451297044490,
                                  -0.887998041597,
                                  0.234528941479,
                                  -0.023185843322};
  const parameter CFX1 = 0.667;
  const parameter CFX2 = 0.8;
  const parameter CFDX1 = 1.24;

  const parameter D_DAMP2 = 0.361;
  const parameter DX_DAMP4_P = 0.232;
  const parameter DX_DAMP4_A = 0.232;
  const parameter B1 = 0.156632;
  const parameter B2 = 0.12083;
  const parameter B3 = 0.5;
  const parameter B4 = MU * MU / K1 - 0.112654;

  const parameter ALPHA_GE = 20.0 / 27.0 + ETA * 5.0 / 3.0;

  // Interpolation function
  num oma = 1.0 - alpha;
  num ief = 0.0;
  if (IINTERP == 0) {
    // SCAN
    if (alpha < 1.0) {
      ief = exp(-CFX1 * alpha / oma);
    } else {
      ief = -CFDX1 * exp(CFX2 / oma);
    }
  } else if (IINTERP == 1) {
    // rSCAN
    if (alpha < 1.0e-13) {
      ief = exp(-CFX1 * alpha / oma);
    } else if (alpha < 2.5) {
      for (int i = 0; i < 8; i++) {
        ief += IE_PARAMS[i] * pow(alpha, i);
      }
    } else {
      ief = -CFDX1 * exp(CFX2 / oma);
    }
  } else {
    printf("ERROR: Unknown IINTERP %d\n", IINTERP);
  }

  // Single orbital enhancement
  num h0x = 1.0 + K0;

  // Slowly varying enhancement
  num del_f2 = 0.0;
  num C2 = 0.0;
  num h1x = 0.0;

  if (IDELFX == 0) {
    // 2nd and 4th order gradient expansion corrections for SCAN interpolation
    num wfac = B4 * pow2(p) * exp(-B4 * p / MU);
    num vfac = B1 * p + B2 * oma * exp(-B3 * pow2(oma));
    num yfac = MU * p + wfac + pow2(vfac);

    h1x = 1.0 + K1 - K1 / (1.0 + yfac / K1);
  } else if (IDELFX == 1 || IDELFX == 2) {
    // 2nd order GE corrections for rSCAN interpolation
    for (int i = 1; i < 8; i++) {
      del_f2 += i * IE_PARAMS[i];
    }
    C2 = -del_f2 * (1.0 - h0x);

    num damp = exp(-pow2(p) / pow(D_DAMP2, 4));
    h1x = 1.0 + K1 - K1 / (1.0 + p * (MU + ALPHA_GE * C2 * damp) / K1);
  } else {
    printf("ERROR: Unknown IDELFX %d\n", IDELFX);
  }

  // Scaling correction
  num gx = 1.0 - exp(-A1 / pow(p, 1.0 / 4.0));

  // 4th order gradient enhancement
  num del_fx = 0.0;
  if (IDELFX == 2) {
    // 4th order correction for rSCAN interpolation (r4SCAN)
    num eta_term = ETA * 3.0 / 4.0 + 2.0 / 3.0;

    num del_f4 = 0.0;
    for (int i = 1; i < 8; i++) {
      del_f4 += i * (i - 1) * IE_PARAMS[i];
    }

    num C_aa = 73.0 / 5000.0 - 0.5 * del_f4 * (h0x - 1.0);
    num C_pa = 511.0 / 13500 - 73.0 / 1500.0 * ETA - del_f2 * (ALPHA_GE * C2 + MU);
    num C_pp = 146.0 / 2025.0 * pow2(eta_term) - 73.0 / 405.0 * eta_term +
               pow2(ALPHA_GE * C2 + MU) / K1;

    num order_1 = C2 * (oma - ALPHA_GE * p);
    num t1 = order_1 + C_aa * pow2(oma) + C_pa * p * oma + C_pp * pow2(p);

    num damp_4_t1 = 2.0 * pow2(alpha) / (1.0 + pow(alpha, 4));
    num damp_4_t2 =
        exp(-pow2(oma) / pow2(DX_DAMP4_A) - pow2(p) / pow(DX_DAMP4_P, 4));
    num damp_4 = damp_4_t1 * damp_4_t2;

    del_fx = t1 * damp_4;
  }

  num fx = (h1x + ief * (h0x - h1x) + del_fx) * gx;

  return fx;
}

template <class num>
static num SCAN_C(const densvars<num> & d,
                  const int IALPHA,
                  const int IINTERP,
                  const int IDELEC) {

  const parameter CFC1 = 0.64;
  const parameter CFC2 = 1.5;
  const parameter CFDC1 = 0.7;
  const parameter IE_PARAMS[8] = {1.0,
                                  -0.64,
                                  -0.4352,
                                  -1.535685604549,
                                  3.061560252175,
                                  -1.915710236206,
                                  0.516884468372,
                                  -0.051848879792};

  const parameter ETA = 1.0e-3;
  const parameter TAU_R = 1.0e-4;
  const parameter A_REG = 1.0e-3;

  const parameter B1C = 0.0285764;
  const parameter B2C = 0.0889;
  const parameter B3C = 0.125541;

  num rs = pow(4.0 * PI * d.n / 3.0, -(1.0 / 3.0));
  num sqrtrs = 0.0;
  if (abs(rs) > 1.0e-16) {
    sqrtrs = sqrt(rs);
  }

  num ds_z = ufunc(d.zeta, 5.0 / 3.0) / 2.0;

  num s = sqrt(d.gnn) / (2.0 * pow(3.0 * PI2, 1.0 / 3.0) * pow(d.n, 4.0 / 3.0));

  num tueg_con = 3.0 / 10.0 * pow(3.0 * PI2, 2.0 / 3.0);
  num tueg = 0.0;
  if (IALPHA == 1) {
    tueg = (tueg_con * pow(d.n, 5.0 / 3.0) + TAU_R) * ds_z;
  } else {
    tueg = tueg_con * pow(d.n, 5.0 / 3.0) * ds_z;
  }

  num tauw = d.gnn / (8.0 * d.n);

  num alpha = 0.0;
  if (IALPHA == 0) {
    // alpha (SCAN)
    if (abs(d.tau - tauw) > 1.0e-14) {
      alpha = (d.tau - tauw) / tueg;
    }

  } else if (IALPHA == 1) {
    // alpha' (rSCAN)
    alpha = (d.tau - tauw) / tueg;
    alpha = pow(alpha, 3) / (pow2(alpha) + A_REG);

  } else if (IALPHA == 2) {
    // \bar{alpha} (r2SCAN, r4SCAN)
    if (abs(d.tau - tauw) > 1.0e-14) {
      alpha = (d.tau - tauw) / (tueg + ETA * tauw);
    }

  } else {
    // Unknown IALPHA choice
    printf("ERROR: Unknown IALPHA %d\n", IALPHA);
  }

  // Interpolation function
  num oma = 1.0 - alpha;
  num ief = 0.0;
  if (IINTERP == 0) {
    // SCAN
    if (alpha < 1.0) {
      ief = exp(-CFC1 * alpha / oma);
    } else {
      ief = -CFDC1 * exp(CFC2 / oma);
    }
  } else if (IINTERP == 1) {
    // rSCAN
    if (alpha < 1.0e-13) {
      ief = exp(-CFC1 * alpha / oma);
    } else if (alpha < 2.5) {
      for (int i = 0; i < 8; i++) {
        ief += IE_PARAMS[i] * pow(alpha, i);
      }
    } else {
      ief = -CFDC1 * exp(CFC2 / oma);
    }
  } else {
    printf("ERROR: Unknown IINTERP %d\n", IINTERP);
  }

  num ec0 = scan_ec0(rs, s, d.zeta, B1C, B2C, B3C);
  num ec1 = scan_ec1(rs, s, d.zeta, IE_PARAMS, ETA, B1C, B2C, B3C, IDELEC);

  num eps_c = (ec1 + ief * (ec0 - ec1)) * d.n;

  return eps_c;
}

template <class num>
static num scan_ec0(const num rs,
                    const num s,
                    const num zeta,
                    const parameter B1C,
                    const parameter B2C,
                    const parameter B3C) {
  const parameter CHI_LD = 0.12802585262625815;

  num eclda = lda_0(rs, B1C, B2C, B3C);

  num dx_z = ufunc(zeta, 4.0 / 3.0) / 2.0;
  num gc_z = (1.0 - 2.363 * (dx_z - 1.0)) * (1 - pow(zeta, 12));

  num w0 = exp(-eclda / B1C) - 1.0;

  num ginf = 1.0 / pow(1.0 + 4.0 * CHI_LD * s * s, 1.0 / 4.0);

  num h0 = B1C * log(1.0 + w0 * (1.0 - ginf));

  return (eclda + h0) * gc_z;
}

template <class num>
static num lda_0(const num rs,
                 const parameter B1C,
                 const parameter B2C,
                 const parameter B3C) {

  return -B1C / (1.0 + B2C * sqrt(rs) + B3C * rs);
}

template <class num>
static num scan_ec1(const num rs,
                    const num s,
                    const num zeta,
                    const parameter IE_PARAMS[8],
                    const parameter ETA,
                    const parameter B1C,
                    const parameter B2C,
                    const parameter B3C,
                    const int IDELEC) {
  const parameter BETA_MB = 0.066725;
  const parameter AFACTOR = 0.1;
  const parameter BFACTOR = 0.1778;
  const parameter GAMMA = 0.031090690869655;
  const parameter AFIX_T = sqrt(PI / 4.0) * pow(9.0 * PI / 4.0, 1.0 / 6.0);
  const parameter D_DAMP2 = 0.361;

  num dx_z = ufunc(zeta, 4.0 / 3.0) / 2.0;
  num gc_z = (1.0 - 2.363 * (dx_z - 1.0)) * (1 - pow(zeta, 12));
  num phi = ufunc(zeta, 2.0 / 3.0) / 2.0;
  num phi3 = pow(phi, 3);

  num sqrtrs = sqrt(rs);

  num eclda0 = lda_0(rs, B1C, B2C, B3C);

  num eclsda1 = 0.0;
  num d_eclsda1_drs = 0.0;
  get_lsda1(rs, sqrtrs, zeta, eclsda1, d_eclsda1_drs);

  num t = AFIX_T * s / (sqrtrs * phi);

  num w1 = exp(-eclsda1 / (GAMMA * phi3)) - 1.0;

  num beta = BETA_MB * (1.0 + AFACTOR * rs) / (1.0 + BFACTOR * rs);
  num y = beta / (GAMMA * w1) * pow2(t);

  num del_y = 0.0;
  if (IDELEC == 0) {
    // No correction for 2nd order GE in correlation (SCAN)
    del_y = 0.0;
  } else if (IDELEC == 1 || IDELEC == 2) {
    // correcting terms for 2nd order GE with rSCAN interpolation
    // Note that IDELEC = 1 is identical to 2

    num p = s * s;
    num ds_z = ufunc(zeta, 5.0 / 3.0) / 2.0;

    num del_f2 = 0.0;
    for (int i = 1; i < 8; i++) {
      del_f2 += i * IE_PARAMS[i];
    }

    num eclsda0 = eclda0 * gc_z;
    num d_eclsda0_drs = gc_z * (B3C + B2C / (2.0 * sqrtrs)) * pow2(eclda0) / B1C;

    num t1 = del_f2 / (27.0 * GAMMA * ds_z * phi3 * w1);
    num t2 = 20.0 * rs * (d_eclsda0_drs - d_eclsda1_drs);
    num t3 = 45.0 * ETA * (eclsda0 - eclsda1);

    num k = t1 * (t2 - t3);

    num damp = exp(-pow2(p) / pow(D_DAMP2, 4));

    del_y = k * p * damp;
  } else {
    printf("ERROR: Unrecognised IDELEC %d\n", IDELEC);
  }

  num g_y = 1.0 / pow(1.0 + 4.0 * (y - del_y), 1.0 / 4.0);

  num h1 = GAMMA * phi3 * log(1.0 + w1 * (1.0 - g_y));

  return eclsda1 + h1;
}

template <class num>
static void get_lsda1(const num rs,
                      const num sqrtrs,
                      const num zeta,
                      num & eclda1,
                      num & d_eclda1_drs) {
  const parameter GAM = 0.51984209978974632953442121455650;
  const parameter FZZ = 8.0 / (9.0 * GAM);
  const parameter p_eu[6] = {
      0.03109070, 0.213700, 7.59570, 3.58760, 1.63820, 0.492940};
  const parameter p_ep[6] = {
      0.015545350, 0.205480, 14.11890, 6.19770, 3.36620, 0.625170};
  const parameter p_alfm[6] = {
      0.01688690, 0.111250, 10.3570, 3.62310, 0.880260, 0.496710};

  num eu = 0.0;
  num deudrs = 0.0;
  gcor2(p_eu, rs, sqrtrs, eu, deudrs);
  num ep = 0.0;
  num depdrs = 0.0;
  gcor2(p_ep, rs, sqrtrs, ep, depdrs);
  num alfm = 0.0;
  num dalfmdrs = 0.0;
  gcor2(p_alfm, rs, sqrtrs, alfm, dalfmdrs);

  num z3 = pow(zeta, 3);
  num z4 = zeta * z3;

  num f = (ufunc(zeta, 4.0 / 3.0) - 2.0) / GAM;

  eclda1 = eu * (1.0 - f * z4) + ep * f * z4 - alfm * f * (1.0 - z4) / FZZ;
  d_eclda1_drs =
      (1.0 - z4 * f) * deudrs + z4 * f * depdrs - (1.0 - z4) * f * dalfmdrs / FZZ;
  return;
}

template <class num>
static void gcor2(const parameter P[6],
                  const num rs,
                  const num sqrtrs,
                  num & GG,
                  num & GGRS) {
  enum IDX { A, A1, B1, B2, B3, B4 };

  num Q0 = -2.0 * P[A] * (1.0 + P[A1] * rs);
  num Q0RS = -2.0 * P[A] * P[A1];

  num Q1 = 2.0 * P[A] * sqrtrs *
           (P[B1] + sqrtrs * (P[B2] + sqrtrs * (P[B3] + P[B4] * sqrtrs)));
  num Q1RS = P[A] * (2.0 * P[B2] + P[B1] / sqrtrs + 3.0 * P[B3] * sqrtrs +
                     4.0 * P[B4] * rs);

  num Q2 = log(1.0 + 1.0 / Q1);
  num Q2RS = -Q1RS / ((1.0 + 1.0 / Q1) * pow2(Q1));

  GG = Q0 * Q2;
  GGRS = Q0 * Q2RS + Q2 * Q0RS;
  return;
}
} // namespace SCAN_eps
