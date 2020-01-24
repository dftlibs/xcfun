.. _using:

Using XCFun
===========

To use the library you should first create one ``xc_functional`` object for each
functional and each thread you want to use. Each functional object can only be
used by one thread at a time. After creating these objects using
:c:func:`xc_new_functional` you should set them up by defining the variables to
differentiate with respect to, as well as the parameters of the functional you
want to use. After having done so you can use :c:func:`xc_eval` to evaluate the
exchange correlation energy and its derivatives.

Example C program that evaluates BLYP to order 2 using alpha/beta variables::

  #include "xcfun.h"

  int main(int argc, char *argv[])
  {
    int derivative_order = 2;
    int nr_points = 1;
    double density[5] = {1,2,3,5,4}; /* na nb ga.ga gb.gb ga.gb */
    double output[21]; /* We have 21 output numbers for derivatives up to order 2 */
    xc_functional fun = xc_new_functional();
    xc_set_mode(fun, XC_VARS_AB);
    xc_set_param(fun, XC_LYPC, 1.0);
    xc_set_param(fun, XC_BECKEX, 1.0);
    xc_eval(fun, derivative_order, nr_points, density, output);
  }

Please see the complete :ref:`api` documentation.

.. _input:

Input, output and units
=======================

The library uses atomic units for all input and output variables.

The XC energy density and derivatives can be evaluated using local spin-up
:math:`(\alpha)` and spin-down :math:`(\beta)` quantities.
In the most general case these are:

    * :math:`n_\alpha` The spin-up electron number density.
    * :math:`n_\beta` The spin-down density.
    * :math:`\sigma_{\alpha \alpha} = \nabla n_\alpha.\nabla n_\alpha` The square magnitude of the spin-up density gradient.
    * :math:`\sigma_{\alpha \beta} = \nabla n_\alpha.\nabla n_\beta` The dot product between the spin-up and spin-down gradient vectors.
    * :math:`\sigma_{\beta \beta} = \nabla n_\beta.\nabla n_\beta` The square magnitude of the spin-down density gradient.
    * :math:`\tau_\alpha = \frac{1}{2} \sum_i |\psi_{i \alpha}|^2` The spin-up Kohn-Sham kinetic energy density.
    * :math:`\tau_\beta` The spin-down Kohn-Sham kinetic energy density. 

Alternatively you can use total density :math:`(n = n_\alpha + n_\beta)` and
spin density :math:`(s = n_\alpha - n_\beta)` variables. These also have
corresponding gradient and kinetic energy components. See :c:func:`xc_set_mode`
below for more information.

The output is given in `graded reverse lexicographical order
<https://en.wikipedia.org/wiki/Monomial_order#Graded_reverse_lexicographic_order>`_.
For example a spin-polarized second order GGA functional will give 21 output elements, starting with the XC energy density. Symbolically we may write this as a list starting with the energy E, followed by five gradient elements
:math:`E_{\alpha} E_{\beta} E_{\sigma_{\alpha \alpha}} E_{\sigma_{\alpha \beta}} E_{\sigma_{\beta \beta}}` 
and 15 second derivatives 
:math:`E_{\alpha \alpha} E_{\alpha \beta} E_{\alpha \sigma_{\alpha \alpha}} ... E_{\beta \beta} E_{\beta \sigma_{\alpha \alpha}} ... E_{\sigma_{\beta \beta} \sigma_{\beta \beta}}` . 
See the function :c:func:`xc_derivative_index` for information about addressing these elements. 



