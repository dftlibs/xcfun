#!/usr/bin/env python

import sys
import numpy
import xcfun


def eval_xcfun_n(fun, data, order=0):

    xcfun.xc_eval_setup(fun, xcfun.XC_N, xcfun.XC_POTENTIAL, order)
    rho = numpy.zeros((data.shape[0], 2))
    rho[:, 0] = data[:, 0]
    rho[:, 1] = 0.0
    out = xcfun.xc_eval(fun, rho)

    return out[:, order]


def main():

    version = xcfun.xcfun_version()
    splash = xcfun.xcfun_splash()

    print("\n Sample use of python interface to xcfun")

    print("\n\nXCFun version ", version)
    print(splash)

    fun_xc_name = 'LDA'
    fun_xc_weigth = 1.0
    fun_xc = xcfun.Functional({fun_xc_name: fun_xc_weigth})

    n_gn_gnn = numpy.zeros((1, 10))
    n_gn_gnn[0, :] = [1.0, 0.2, 0.2, 0.2, 0, 0, 0, 0, 0, 0]

    n = numpy.zeros((1, 1))
    n[0, :] = [1.0]

    result = fun_xc.eval_potential_n(n)

    energy = result[0, 0]
    potential = result[0, 1]

    print("Functional used: ", fun_xc_name, " weight ", fun_xc_weigth)
    print("        density: ", n)
    print("         energy: ", energy)
    print("      potential: ", potential)

    resut_energy = -0.8101513
    resut_potential = -1.06468341

    tol = 1.0e-7
    energy_diff = energy - resut_energy
    potential_diff = potential - resut_potential

    if ((abs(energy_diff) > tol) or (abs(potential_diff) > tol)):
        sys.exit(-1)


if __name__ == '__main__':
    main()
