# -*- coding: utf-8 -*-
#
# XCFun, an arbitrary order exchange-correlation library
# Copyright (C) 2019 Ulf Ekström and contributors.
#
# This file is part of XCFun.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# For information on the complete list of contributors to the
# XCFun library, see: <https://xcfun.readthedocs.io/>
#

import numpy

from .xcfun_swig import *


class XCFunException(Exception):
    def __init__(self, msg, code=1):
        Exception.__init__(self, msg)
        self.error_code = code


def xc_new_functional():
    return xc_new_functional_not_macro(XCFUN_API_VERSION)


def xc_eval_setup(fun, var, mode, order):
    code = xc_eval_setup_swig(fun, var, mode, order)
    if code > 0:
        raise XCFunException("Invalid options in xc_eval_setup", code)


def xc_eval(fun, density):

    dens_len = xc_input_length(fun)
    output_len = xc_output_length(fun)

    if not ((density.shape[-1] == dens_len)):
        raise XCFunException('wrong dimension of density argument')

    if len(density.shape) == 1:
        output = numpy.zeros(output_len)
        xc_eval_swig(fun, density, output)

    elif len(density.shape) == 2:
        nr_points = density.shape[0]
        output = numpy.zeros((nr_points, output_len))
        xc_eval_vec_swig(fun, density, output)

    else:
        raise XCFunException('wrong shape of density argument')

    return output


class Functional(object):
    """
    Python class representing an XCFun functional.
    """

    def __init__(self, funcdict):
        """
        Initialize the functional.

        funclist is a dictonary containing functional names and corresponding coefficients.
        Examples:
            lda_func = Functional({'LDA': 1.0})
            blyp_func = Functional({'BeckeX': 1.0, 'LYPC': 1.0})
        """
        self._func = xc_new_functional()
        for name, weight in funcdict.items():
            ret = xc_set(self._func, name, weight)
            if not ret == 0:
                raise XCFunException('unknown functional selected')

    def __del__(self):
        xc_free_functional(self._func)

    @property
    def type(self):
        if xc_is_metagga(self._func):
            return 2
        elif xc_is_gga(self._func):
            return 1
        else:
            return 0

    def eval_potential_n(self, density, densgrad=None, denshess=None):
        """
        Evaluate the xc potential (spin-compensated case).

        input:
            density: 1D-numpy.array[0:nr_of_points],
                density at grid points
            densgrad: 2D-numpy.array[0:nr_of_points, 0:3],
                density gradient at grid points (1st index),
                2nd index: 0 - x, 1 - y, 2 - z
                - only required for GGA functionals
            denshess: 2D-numpy.array[0:nr_of_points, 0:6]
                density Hessian at grid points (1st index),
                2nd index: 0 - xx, 1 - xy, 2 - xz, 3 - yy, 4 - yz, 5 - zz
                - only required for GGA functionals
        output:
            return value: 2D-numpy array[0:nr_of_points,0:2]
            2nd index: 0 - energy density, 1 - total xc potential
        """

        if xc_is_metagga(self._func):
            raise XCFunException('xc potential not supported for meta-GGAs')
        elif xc_is_gga(self._func):
            if (densgrad is None) or (denshess is None):
                raise XCFunException('Density gradient and Laplacian required for GGA potential')

            xc_eval_setup(self._func, XC_N_2ND_TAYLOR, XC_POTENTIAL, 1)

            if not (len(density.shape) == 1):
                raise XCFunException('Wrong shape of density argument in eval_potential_n')
            nr_points = density.size

            if not (densgrad.shape == (nr_points, 3)):
                raise XCFunException('Wrong shape of densgrad argument in eval_potential_n')

            if not (denshess.shape == (nr_points, 6)):
                raise XCFunException('Wrong shape of denshess argument in eval_potential_n')

            dens = numpy.zeros((density.size, 10))
            dens[:, 0] = density[:]
            dens[:, 1:4] = densgrad[:, 0:3]
            dens[:, 4:10] = denshess[:, 0:6]

        else:
            xc_eval_setup(self._func, XC_N, XC_POTENTIAL, 1)

            dens = density.reshape((density.size, 1))

        return xc_eval(self._func, dens)

    def eval_potential_ab(self, density, densgrad=None, denshess=None):
        """
        Evaluate the xc potential (spin-resolved case).

        input:
            density: 1D-numpy.array[0:nr_of_points,0:2],
                density at grid points (1st index),
                2nd index: 0 - alpha density, 1 - beta density
            densgrad: 2D-numpy.array[0:nr_of_points, 0:3, 0:2],
                density gradient at grid points (1st index),
                2nd index: 0 - x, 1 - y, 2 - z,
                3rd index: 0 - alpha density gradien, 1 - beta density gradient,
                - only required for GGA functionals
            denshess: 2D-numpy.array[0:nr_of_points, 0:6]
                density Hessian at grid points (1st index),
                2nd index: 0 - xx, 1 - xy, 2 - xz, 3 - yy, 4 - yz, 5 - zz,
                3rd index, 0 - alpha density Hessian, 1 - beta density Hessian
                - only required for GGA functionals
        output:
            return value: 2D-numpy array[0:nr_of_points,0:3]
            2nd index: 0 - energy density, 1 - alpha xc potential, 2 - beta xc potential
        """

        if xc_is_metagga(self._func):
            raise XCFunException('xc potential not supported for meta-GGAs')
        elif xc_is_gga(self._func):
            if (densgrad is None) or (denshess is None):
                raise XCFunException('Density gradient and Laplacian required for GGA potential')

            xc_eval_setup(self._func, XC_A_B_2ND_TAYLOR, XC_POTENTIAL, 1)

            if not (len(density.shape) == 2):
                raise XCFunException('Wrong shape of density argument in eval_potential_ab')
            nr_points = density.shape[0]

            if not (density.shape == (nr_points, 2)):
                raise XCFunException('Wrong shape of density argument in eval_potential_n')

            if not (densgrad.shape == (nr_points, 3, 2)):
                raise XCFunException('Wrong shape of densgrad argument in eval_potential_n')

            if not (denshess.shape == (nr_points, 6, 2)):
                raise XCFunException('Wrong shape of denshess argument in eval_potential_n')

            dens = numpy.zeros((density.shape[0], 20))
            dens[:, 0] = density[:, 0]
            dens[:, 1:4] = densgrad[:, 0:3, 0]
            dens[:, 4:10] = denshess[:, 0:6, 0]
            dens[:, 10] = density[:, 1]
            dens[:, 11:14] = densgrad[:, 0:3, 1]
            dens[:, 14:20] = denshess[:, 0:6, 1]

        else:
            xc_eval_setup(self._func, XC_A_B, XC_POTENTIAL, 1)

            dens = density

        return xc_eval(self._func, dens)
