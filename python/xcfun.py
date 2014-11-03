#!/usr/bin/env python

import numpy
from xcfun_swig import *

class XcFunException(Exception):
    def __init__ (self, msg, code) :
       Exception.__init__(self, msg)
       self.error_code = code

def xc_new_functional() :
    return xc_new_functional_not_macro(XCFUN_API_VERSION)

def xc_eval_setup(fun, var, mode, order) :
    code = xc_eval_setup_swig (fun, var, mode, order)
    if code > 0 :
        raise XcFunException("Invalid options in xc_eval_setup", code)

def xc_eval(fun, density) :

     #dens_len = xc_input_length(fun)
     output_len = xc_output_length(fun)

     #if not (density.shape[-1] == dens_len) :
     #    raise XcFunException('wrong dimension of density argument')

     if len(density.shape) == 1 :
         output = numpy.zeros(output_len)
         xc_eval_swig (fun,  density, output)

     elif len(density.shape) == 2 :
         nr_points = density.shape[0]
         output = numpy.zeros((nr_points, output_len))
         xc_eval_vec_swig (fun, density, output)

     else :
         raise XcFunException('wrong shape of density argument')

     return output


