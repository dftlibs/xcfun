.. _using:

Using XCFun
===========

To use the library, you will need to:

- Link your executable to it. Either using the static, ``libxcfun.a`` or shared,
  ``libxcfun.so``, version.
- For C/C++ hosts, include the header file ``xcfun.h`` where appropriate:

.. code-block:: cpp

   #include "XCFun/xcfun.h"

- For Fortran hosts, compile the ``xcfun.f90`` source file together with your
  sources. This will allow using the Fortran/C interoperability layer with:

.. code-block:: fortran

   use xcfun


Installing using Spack
----------------------

XCFun can be installed in a `Spack
<https://www.spack.io/>`_ environment with::

  spack env create myenv
  spack env activate myenv
  spack install xcfun


.. _integration:

Integration with your build system
----------------------------------

The set up of the build system for you code will change the details on how to
achieve the points above. In the following, we provide minimalistic instructions
for codes that use either `CMake <https://cmake.org/>`_ as their build system
generator or plain ``Makefile``.

.. _cmake-integration:

CMake as build system
~~~~~~~~~~~~~~~~~~~~~

.. note:: You can find **complete**, **standalone** examples for C, C++, and
          Fortran in the ``examples`` folder.

If you use CMake as your build system, adding the command:

.. code-block:: cmake

   find_package(XCFun CONFIG)

in your ``CMakeLists.txt`` will let CMake search for an XCFun installation. CMake will honor the hint variable:

.. code-block:: bash

   -DXCFun_DIR=<install-prefix>/share/cmake/XCFun

and set up the target ``XCFun::xcfun`` for you to link your target against:

.. code-block:: cmake

   target_link_libraries(<your-target-name>
     PRIVATE
       XCFun::xcfun
     )

For Fortran hosts the ``xcfun.f90`` will have to be compiled too. The following addition suffices:

.. code-block:: cmake

   target_sources(<your-target-name>
     PRIVATE
       ${XCFun_Fortran_SOURCES}
     )

.. _other-integration:

Other build systems
~~~~~~~~~~~~~~~~~~~

You will need to set:

- The linker path:

.. code-block:: bash

   -L<install-prefix>/lib64 -lxcfun

note that on some systems it might be ``lib`` rather than ``lib64``.

- For C/C++ codes, the include path:

.. code-block:: bash

   -I<install-prefix>/include

- For Fortran codes, the location of the Fortran/C interoperability source file ``xcfun.f90``:

.. code-block:: bash

   <install-prefix>/include/XCFun/xcfun.f90

.. _interfacing:

Writing an interface 
--------------------

.. note:: Please, read the full :ref:`api` documentation for a complete overview. 

The library exposes an opaque type, :cpp:type:`xcfun_t`, through which you can
obtain the exchange-correlation functional derivatives to the desired order.
To do so:

1. Create one :cpp:type:`xcfun_t` object. There should be **only one** such
   object per thread and per XC functional. In C/C++ this is achieved with:

   .. code-block:: c

      xcfun_t * fun = xcfun_new();

   whereas in Fortran:

   .. code-block:: fortran

      use, intrinsic :: iso_c_binding

      type(c_ptr) :: fun

      fun = xcfun_new()
 
2. The :cpp:type:`xcfun_t` object is now a blank slate. You will need to set the
   exchange-correlation admixture, *i.e.* which functional and which amount to
   use for exchange and correlation. This is achieved with calls to
   :cpp:func:`xcfun_set`:

   .. code-block:: c

      int ierr = 0;
      ierr = xcfun_set(fun, "blyp", 0.9);
      ierr = xcfun_set(fun, "pbec", 0.1);

   We have now set up the BLYP GGA functional.

3. Next, you will have to set up the evaluation strategy, *i.e.* which variables
   will be passed in as input to the functional, which outputs are expected, and
   the order of the derivatives to return upon evaluation. This can be done by
   calling :cpp:func:`xcfun_eval_setup`:

   .. code-block:: c

      ierr = xcfun_eval_setup(fun, XC_A_B_AX_AY_AX_BX_BY_BZ, XC_PARTIAL_DERIVATIVES, 1);

   The convenience function :cpp:func:`xcfun_user_eval_setup` is also available.
   With this set up, we will obtain functional derivatives of the BLYP
   functional up to first order, using :math:`\alpha` and :math:`\beta`
   variables and partial derivatives.

4. We are now ready to run the evaluation and for this you will have to allocate
   a properly sized chunk of memory. The function
   :cpp:func:`xcfun_output_length` will return how large such a scratch array
   has to be:

   .. code-block:: c

      int nout = xcfun_output_length(fun); 

      double * output = malloc(sizeof(double) * nout);

5. Finally, we proceed to the evaluation. We call :cpp:func:`xcfun_eval` with an
   array of density values:

   .. code-block:: c

      xcfun_eval(fun, d_elements, output);

6. The important last step is to clean up the used heap memory.
   :cpp:func:`xcfun_delete` is the function to call:

   .. code-block:: c

      free(output);
      xcfun_delete(fun);

.. _input:

Input, output and units
~~~~~~~~~~~~~~~~~~~~~~~

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
corresponding gradient and kinetic energy components. See :cpp:func:`xcfun_set`
below for more information.

The output is given in `graded reverse lexicographical order
<https://en.wikipedia.org/wiki/Monomial_order#Graded_reverse_lexicographic_order>`_.
For example a spin-polarized second order GGA functional will give 21 output elements, starting with the XC energy density. Symbolically we may write this as a list starting with the energy E, followed by five gradient elements
:math:`E_{\alpha} E_{\beta} E_{\sigma_{\alpha \alpha}} E_{\sigma_{\alpha \beta}} E_{\sigma_{\beta \beta}}` 
and 15 second derivatives 
:math:`E_{\alpha \alpha} E_{\alpha \beta} E_{\alpha \sigma_{\alpha \alpha}} ... E_{\beta \beta} E_{\beta \sigma_{\alpha \alpha}} ... E_{\sigma_{\beta \beta} \sigma_{\beta \beta}}` . 
