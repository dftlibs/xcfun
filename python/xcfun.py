#!/usr/bin/env python

import numpy
from xcfun_swig import *

class XcFunException(Exception):
    def __init__ (self, msg, code=1) :
       Exception.__init__(self, msg)
       self.error_code = code

def xc_new_functional() :
    return xc_new_functional_not_macro(XCFUN_API_VERSION)

def xc_eval_setup(fun, var, mode, order) :
    code = xc_eval_setup_swig (fun, var, mode, order)
    if code > 0 :
        raise XcFunException("Invalid options in xc_eval_setup", code)

def xc_eval(fun, density) :

     dens_len = xc_input_length(fun)
     output_len = xc_output_length(fun)

     if not ((density.shape[-1] == dens_len)) :
         raise XcFunException('wrong dimension of density argument')

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

class Functional(object) :
    """
    Python class representing an XCFun functional.
    """

    def __init__ (self, funcdict) :
        """
        Initialize the functional.

        funclist is a dictonary containing functional names and corresponding coefficients.
        Examples:
            lda_func = Functional({'LDA': 1.0})
            blyp_func = Functional({'BeckeX': 1.0, 'LYPC': 1.0})
        """
        self._func = xc_new_functional()
        for (name, weight) in funcdict.iteritems():
            xc_set(self._func, name, weight)

    def __del__ (self) :
        xc_free_functional(self._func)

    def eval_potential_n (self, density, densgrad=None, denslapl=None) :
        """
        Evaluate the xc potential (spin-compensated case).

        input:
        output:
        """

        if xc_is_metagga(self._func) :
            raise XcFunException('xc potential not supported for meta-GGAs')
        elif xc_is_gga(self._func) :
            if (densgrad is None) or (denslapl is None) :
                raise XcFunException('Density gradient and Laplacian required for GGA xc potential')

            xc_eval_setup(self._func, XC_N_2ND_TAYLOR, XC_POTENTIAL, 1)

            dens = numpy.zeros((density.size, 10))
            dens[:,0] = density[:]
            dens[:,1:4] = densgrad[:,0:3]
            dens[:,4:10] = denslapl[:,0:6]

        else :
            xc_eval_setup(self._func, XC_N, XC_POTENTIAL, 1)

            dens = density.reshape((density.size, 1))
 
        return xc_eval(self._func, dens)

    def eval_potential_ab (self, density, densgrad=None, denslapl=None) :
        """
        Evaluate the xc potential (spin-resolved case).

        input:
        output:
        """

        if xc_is_metagga(self._func) :
            raise XcFunException('xc potential not supported for meta-GGAs')
        elif xc_is_gga(self._func) :
            if (densgrad is None) or (denslapl is None) :
                raise XcFunException('Density gradient and Laplacian required for GGA xc potential')

            xc_eval_setup(self._func, XC_A_B_2ND_TAYLOR, XC_POTENTIAL, 1)

            dens = numpy.zeros((density.shape[0], 20))
            dens[:,0] = density[:,0]
            dens[:,1:4] = densgrad[:,0:3,0]
            dens[:,4:10] = denslapl[:,0:6,0]
            dens[:,10] = density[:,1]
            dens[:,11:14] = densgrad[:,0:3,1]
            dens[:,14:20] = denslapl[:,0:6,1]

        else :
            xc_eval_setup(self._func, XC_A_B, XC_POTENTIAL, 1)

            dens = density

        return xc_eval(self._func, dens)

