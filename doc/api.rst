.. _api:

***
API
***

The library is written in C++, but can also be directly used in a C project or Fortran project. The C interface is described in ``include/xcfun.h``, while the Fortran interface is described in the module file ``fortran/xcfun_module.f90``. This documentation tries to describe both the C and Fortran API at the same time, even though small differences may exist between the two interfaces. 

.. _setup:

Setup and testing
=================


.. c:function:: const char *xcfun_splash(void)

Return a multi-line string describing the library. Please print this string so that your users find the right citation for the library. The Fortran version is ``xcfun_splash(text)``, and the message is put in the string text.


.. c:function:: double xcfun_version(void)

Return a double precision version number of the library.


.. c:function:: int xcfun_test(void)

Run all internal tests and return the number of failed tests.

.. _create_destroy:

Creating and destroying functional objects
==========================================

.. c:function:: xc_functional xc_new_functional()


Create a new functional object. The C version returns an object of type ``xc_functional``. The Fortran version returns an integer. The creation of this object may be rather slow; create an object once for each calculation, not once for each grid point.


.. c:function:: void xc_free_functional(xc_functional fun)

Release the memory associated with functional (previously allocated by :c:func:`xc_new_functional`). 

.. _def_func:

Defining a functional
=====================



.. c:function:: void xc_eval(xc_functional fun, const double *density, double *result)

   Evaluate a functional and its derivatives.

   :param fun: Functional object previously created by :c:func:`xc_new_functional` .
   :param density: Input density parameters, in the order previously defined by :c:func:`xc_eval_setup` .
   :param result: Output values, depending the the mode and order.


.. c:function:: int xc_derivative_index(xc_functional fun, const int derivative[])

Given an integer array of "exponents" this function returns the index into the output array of :c:func:`xc_eval` where the corresponding derivative is located. Note that this position depends on the functional type and mode. Alternatively you can, for low orders, use the predefined constants of the form ``XC_D01000`` (which corresponds to the first derivative with respect to the second variable, for setups with five variables). To use these predefined constants you must therefore know the type of your functional beforehand. For LDA in a two variable mode the corresponding constant would be ``XC_D01``, but this is a different index from ``XC_D01000``. 


.. c:function:: int xc_output_length(xc_functional fun)

Return the number of output coefficients computed by :c:func:`xc_eval` if functional is evaluate up to order. Note that all derivatives up to order are calculated, not only those of the particular order. 

.. _func_list:

List of existing functions in XCFun
===================================


* :c:func:`xcfun_version`

* :c:func:`xcfun_splash`

* :c:func:`xcfun_test`

* :c:func:`xc_new_functional`

* :c:func:`xc_free_functional`

* :c:func:`xc_eval`

* :c:func:`xc_derivative_index`

* :c:func:`xc_output_length`

**typedef struct xc_functional_obj * xc_functional**


.. c:function:: const char *xc_enumerate_parameters(int n)


.. c:function:: const char *xc_enumerate_aliases(int n)


.. c:function:: int xc_set(xc_functional fun, const char *name, double value)


.. c:function:: int xc_get(xc_functional fun, const char *name, double *value)


.. c:function:: const char *xc_describe_short(const char *name)


.. c:function:: const char *xc_describe_long(const char *name)


.. c:function:: int xc_eval_setup(xc_functional fun, enum xc_vars vars, enum xc_mode mode, int order)


.. c:function:: void xc_eval_vec(xc_functional fun, int nr_points, const double *density, int density_pitch, double *result, int result_pitch)




