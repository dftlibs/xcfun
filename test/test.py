#!/usr/bin/env python

import re
import math, numpy
import xcfun

def eval_xcfun_n(fun, data, order=0):

    xcfun.xc_eval_setup(fun, xcfun.XC_N, xcfun.XC_POTENTIAL, order)
    rho = numpy.zeros((data.shape[0],2))
    rho[:,0] = data[:,0]
    rho[:,1] = 0.0
    out = xcfun.xc_eval(fun, rho)

    return out[:,order]


def main():

    version = xcfun.xcfun_version()
    splash  = xcfun.xcfun_splash()

    print "\n Sample use of python interface to xcfun"

    print "\n\nXCFun version ",version
    print splash

    fun_xc = xcfun.xc_new_functional()
    fun_xc_name = "lda"
    fun_xc_weigth = 1.0
    xcfun.xc_set(fun_xc, fun_xc_name, fun_xc_weigth);

    n_gn_gnn = numpy.zeros((1,10))
    n_gn_gnn[0,:] = [1.0, 0.2, 0.2, 0.2,  0, 0, 0, 0, 0, 0]

    energy    = eval_xcfun_n(fun_xc, n_gn_gnn, order=0)
    potential = eval_xcfun_n(fun_xc, n_gn_gnn, order=1)

    print "Functional used: ",fun_xc_name," weight ",fun_xc_weigth
    print "        density: ",n_gn_gnn  
    print "         energy: ",energy
    print "      potential: ",potential

main()
