! XCFun, an arbitrary order exchange-correlation library
! Copyright (C) 2020 Ulf Ekström and contributors.
!
! This file is part of XCFun.
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
! For information on the complete list of contributors to the
! XCFun library, see: <https://xcfun.readthedocs.io/>

module xcfun

  use, intrinsic :: iso_c_binding

  implicit none
!
  integer, parameter :: XCFUN_API_VERSION = 2

! Evaluation modes
  integer, parameter :: XC_PARTIAL_DERIVATIVES = 1
  integer, parameter :: XC_POTENTIAL = 2
  integer, parameter :: XC_CONTRACTED = 3

! Variable combinations
  integer, parameter :: XC_A = 0
  integer, parameter :: XC_N = 1
  integer, parameter :: XC_A_B = 2
  integer, parameter :: XC_N_S = 3
  integer, parameter :: XC_A_GAA = 4
  integer, parameter :: XC_N_GNN = 5
  integer, parameter :: XC_A_B_GAA_GAB_GBB = 6
  integer, parameter :: XC_N_S_GNN_GNS_GSS = 7
  integer, parameter :: XC_A_GAA_LAPA = 8
  integer, parameter :: XC_A_GAA_TAUA = 9
  integer, parameter :: XC_N_GNN_LAPN = 10
  integer, parameter :: XC_N_GNN_TAUN = 11
  integer, parameter :: XC_A_B_GAA_GAB_GBB_LAPA_LAPB = 12
  integer, parameter :: XC_A_B_GAA_GAB_GBB_TAUA_TAUB = 13
  integer, parameter :: XC_N_S_GNN_GNS_GSS_LAPN_LAPS = 14
  integer, parameter :: XC_N_S_GNN_GNS_GSS_TAUN_TAUS = 15
  integer, parameter :: XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB = 16
  integer, parameter :: XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB = 17
  integer, parameter :: XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS = 18
  integer, parameter :: XC_A_AX_AY_AZ = 19
  integer, parameter :: XC_A_B_AX_AY_AZ_BX_BY_BZ = 20
  integer, parameter :: XC_N_NX_NY_NZ = 21
  integer, parameter :: XC_N_S_NX_NY_NZ_SX_SY_SZ = 22
  integer, parameter :: XC_A_AX_AY_AZ_TAUA = 23
  integer, parameter :: XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB = 24
  integer, parameter :: XC_N_NX_NY_NZ_TAUN = 25
  integer, parameter :: XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS = 26
  integer, parameter :: XC_A_2ND_TAYLOR = 27
  integer, parameter :: XC_A_B_2ND_TAYLOR = 28
  integer, parameter :: XC_N_2ND_TAYLOR = 29
  integer, parameter :: XC_N_S_2ND_TAYLOR = 30

! Indices into the output array of derivatives. Fortran numbering.
  integer, parameter :: XC_D0 = 1
  integer, parameter :: XC_D1 = 2
  integer, parameter :: XC_D2 = 3

! all lda derivatives up to order 4
  integer, parameter :: XC_D00 = 1
  integer, parameter :: XC_D10 = 2
  integer, parameter :: XC_D01 = 3
  integer, parameter :: XC_D20 = 4
  integer, parameter :: XC_D11 = 5
  integer, parameter :: XC_D02 = 6
  integer, parameter :: XC_D30 = 7
  integer, parameter :: XC_D21 = 8
  integer, parameter :: XC_D12 = 9
  integer, parameter :: XC_D03 = 10
  integer, parameter :: XC_D40 = 11
  integer, parameter :: XC_D31 = 12
  integer, parameter :: XC_D22 = 13
  integer, parameter :: XC_D13 = 14
  integer, parameter :: XC_D04 = 15

! gga up to order 3
! only derivatives that are nonzero at closed-shell reference
! with XC_VARS_NS are here
  integer, parameter :: XC_D00000 = 1
  integer, parameter :: XC_D10000 = 2
  integer, parameter :: XC_D00100 = 4
  integer, parameter :: XC_D00001 = 6
  integer, parameter :: XC_D20000 = 7
  integer, parameter :: XC_D10100 = 9
  integer, parameter :: XC_D10001 = 11
  integer, parameter :: XC_D02000 = 12
  integer, parameter :: XC_D01010 = 14
  integer, parameter :: XC_D00200 = 16
  integer, parameter :: XC_D00101 = 18
  integer, parameter :: XC_D00020 = 19
  integer, parameter :: XC_D30000 = 22
  integer, parameter :: XC_D20100 = 24
  integer, parameter :: XC_D12000 = 27
  integer, parameter :: XC_D11010 = 29
  integer, parameter :: XC_D10200 = 31
  integer, parameter :: XC_D10020 = 34
  integer, parameter :: XC_D02100 = 38
  integer, parameter :: XC_D01110 = 42
  integer, parameter :: XC_D00300 = 47
  integer, parameter :: XC_D00120 = 50

! tau mgga
! so far only linear response
! list under construction ...
  integer, parameter :: XC_D0000000 = 1
  integer, parameter :: XC_D1000000 = 2
  integer, parameter :: XC_D0010000 = 4
  integer, parameter :: XC_D0000100 = 6
  integer, parameter :: XC_D0000010 = 7
  integer, parameter :: XC_D0000001 = 8
  integer, parameter :: XC_D2000000 = 9
  integer, parameter :: XC_D1010000 = 11
  integer, parameter :: XC_D1000010 = 14
  integer, parameter :: XC_D1000001 = 15
  integer, parameter :: XC_D0200000 = 16
  integer, parameter :: XC_D0101000 = 18
  integer, parameter :: XC_D0100010 = 20
  integer, parameter :: XC_D0100001 = 21
  integer, parameter :: XC_D0020000 = 22
  integer, parameter :: XC_D0010010 = 25
  integer, parameter :: XC_D0010001 = 26
  integer, parameter :: XC_D0002000 = 27
  integer, parameter :: XC_D0001010 = 29
  integer, parameter :: XC_D0001001 = 30
  integer, parameter :: XC_D0000110 = 32
  integer, parameter :: XC_D0000101 = 33
  integer, parameter :: XC_D0000020 = 34
  integer, parameter :: XC_D0000011 = 35
  integer, parameter :: XC_D0000002 = 36

  private :: fstring_to_carray
  private :: text_handler

  interface
    function xcfun_version() result(version) &
      bind(C)
      import
      real(c_double) :: version
    end function

    function xcfun_splash_C() result(text) &
      bind(C, name="xcfun_splash")
      import
      type(c_ptr) :: text
    end function

    function xcfun_test() result(err) &
      bind(C)
      import
      integer(c_int) :: err
    end function

    function xcfun_is_compatible_library() result(is_compatible) &
      bind(C)
      import
      logical(c_bool) :: is_compatible
    end function

    function xcfun_enumerate_parameters_C(n) result(text) &
      bind(C, name="xcfun_enumerate_parameters")
      import
      integer(c_int), intent(in), value :: n
      type(c_ptr) :: text
    end function

    function xcfun_enumerate_aliases_C(n) result(text) &
      bind(C, name="xcfun_enumerate_aliases")
      import
      integer(c_int), intent(in), value :: n
      type(c_ptr) :: text
    end function

    function xcfun_describe_short_C(name) result(text) &
      bind(C, name="xcfun_describe_short")
      import
      character(kind=c_char, len=1), intent(in) :: name(*)
      type(c_ptr) :: text
    end function

    function xcfun_describe_long_C(name) result(text) &
      bind(C, name="xcfun_describe_long")
      import
      character(kind=c_char, len=1), intent(in) :: name(*)
      type(c_ptr) :: text
    end function

    function xcfun_new() result(fun) &
      bind(C)
      import
      type(c_ptr) :: fun
    end function

    subroutine xcfun_delete(fun) &
      bind(C)
      import
      type(c_ptr), value :: fun
    end subroutine

    function xcfun_set_C(fun, name, val) result(err) &
      bind(C, name="xcfun_set")
      import
      type(c_ptr), value :: fun
      character(kind=c_char, len=1), intent(in) :: name(*)
      real(c_double), intent(in), value :: val
      integer(c_int) :: err
    end function

    function xcfun_get_C(fun, name, val) result(err) &
      bind(C, name="xcfun_get")
      import
      type(c_ptr), value :: fun
      character(kind=c_char, len=1), intent(in) :: name(*)
      real(c_double), intent(inout) :: val(*)
      integer(c_int) :: err
    end function

    subroutine xcfun_eval_vec_C(fun, nr_points, density, d_pitch, res, r_pitch) &
         bind(C, name="xcfun_eval_vec")
      import
      type(c_ptr), intent(in), value :: fun
      integer(c_int), intent(in), value :: nr_points
      real(c_double), intent(in) :: density(*)
      integer(c_int), intent(in), value :: d_pitch
      real(c_double), intent(inout) :: res(*)
      integer(c_int), intent(in), value :: r_pitch
    end subroutine

    function xcfun_is_gga(fun) result(is_gga) &
      bind(C)
      import
      type(c_ptr), value :: fun
      logical(c_bool) :: is_gga
    end function

    function xcfun_is_metagga(fun) result(is_metagga) &
      bind(C)
      import
      type(c_ptr), value :: fun
      logical(c_bool) :: is_metagga
    end function

    function xcfun_input_length(fun) result(length) &
      bind(C)
      import
      type(c_ptr), intent(in), value :: fun
      integer(c_int) :: length
    end function

    function xcfun_output_length(fun) result(length) &
      bind(C)
      import
      type(c_ptr), intent(in), value :: fun
      integer(c_int) :: length
    end function

    function xcfun_derivative_index(fun, derivative) result(index) &
      bind(C)
      import
      type(c_ptr), intent(in), value :: fun
      integer(c_int), intent(in) :: derivative(*)
      integer(c_int) :: index
    end function
  end interface

  interface xcfun_eval_setup
    function xcfun_eval_setup(fun, vars, mode, order) result(err) &
         bind(C)
      import
      type(c_ptr), value :: fun
      integer(c_int), intent(in), value :: vars
      integer(c_int), intent(in), value :: mode
      integer(c_int), intent(in), value :: order
      integer(c_int) :: err
    end function

    function xcfun_user_eval_setup(fun, order, func_type, dens_type, mode_type, &
         laplacian, kinetic, current, explicit_derivatives) result(err) &
         bind(C)
      import
      type(c_ptr), value :: fun
      integer(c_int), intent(in), value :: order
      integer(c_int), intent(in), value :: func_type
      integer(c_int), intent(in), value :: dens_type
      integer(c_int), intent(in), value :: mode_type
      integer(c_int), intent(in), value :: laplacian
      integer(c_int), intent(in), value :: kinetic
      integer(c_int), intent(in), value :: current
      integer(c_int), intent(in), value :: explicit_derivatives
      integer(c_int) :: err
    end function
  end interface

  interface xcfun_eval
    subroutine xcfun_eval(fun, density, res) &
         bind(C)
      import
      type(c_ptr), intent(in), value :: fun
      real(c_double), intent(in) :: density(*)
      real(c_double), intent(inout) :: res(*)
    end subroutine

    module procedure xcfun_eval_vec
  end interface

contains
  ! \brief Convert a Fortran string into a C string.
  ! \param[in] string_f03 a Fortran character string.
  ! \return array_c Null-terminated C string in a character array.
  pure function fstring_to_carray(string_f03) result(array_c)
    character(len=*), intent(in) :: string_f03
    character(kind=c_char, len=1) :: array_c(len(string_f03)+1)

    integer :: i

    do i = 1, len(string_f03)
       array_c(i) = string_f03(i:i)
    end do
    array_c(i) = c_null_char
  end function

  function text_handler(C_text, length) result(F_text)
    integer(c_int), intent(in) :: length
    type(c_ptr), intent(in) :: C_text
    character(kind=c_char, len=length) :: F_text

    character(kind=c_char), pointer, dimension(:) :: msg_array
    character(kind=c_char, len=length) :: msg
    integer(c_int) :: msg_length
    integer :: i

    call c_f_pointer(C_text, msg_array, [ length ])

    do i = 1, length
       msg(i:i+1) = msg_array(i)
    end do

    msg_length = len_trim(msg(1:index(msg, c_null_char)))

    F_text = msg(1:msg_length-1)
  end function

  function xcfun_splash() result(text)
    character(kind=c_char, len=5000) :: text

    text = text_handler(xcfun_splash_C(), 5000)
  end function

  function xcfun_enumerate_parameters(n) result(text)
    integer(c_int), intent(in) :: n
    character(kind=c_char, len=5000) :: text

    text = text_handler(xcfun_enumerate_parameters_C(n), 5000)
  end function

  function xcfun_enumerate_aliases(n) result(text)
    integer(c_int), intent(in) :: n
    character(kind=c_char, len=10000) :: text

    text = text_handler(xcfun_enumerate_aliases_C(n), 10000)
  end function

  function xcfun_describe_short(name) result(text)
    character(kind=c_char, len=1), intent(in) :: name(*)
    character(kind=c_char, len=10000) :: text

    text = text_handler(xcfun_describe_short_C(name), 10000)
  end function

  function xcfun_describe_long(name) result(text)
    character(kind=c_char, len=1), intent(in) :: name(*)
    character(kind=c_char, len=10000) :: text

    text = text_handler(xcfun_describe_long_C(name), 10000)
  end function

  function xcfun_set(fun, param, val) result(err)
    type(c_ptr), value :: fun
    character(kind=c_char, len=*), intent(in) :: param
    real(c_double), intent(in) :: val
    integer(c_int) :: err

    err = xcfun_set_C(fun, fstring_to_carray(param), val)
  end function

  function xcfun_get(fun, param, val) result(err)
    type(c_ptr), value :: fun
    character(kind=c_char, len=*), intent(in) :: param
    real(c_double), intent(inout) :: val(*)
    integer(c_int) :: err

    err = xcfun_get_C(fun, fstring_to_carray(param), val)
  end function

  subroutine xcfun_eval_vec(fun, nr_points, density, res)
    type(c_ptr), intent(in), value :: fun
    integer(c_int), intent(in) :: nr_points
    real(c_double), intent(in) :: density(:, :)
    real(c_double), intent(inout) :: res(:, :)

    integer(c_int) :: d_pitch
    integer(c_int) :: r_pitch

    d_pitch = size(density, kind=c_int)
    r_pitch = size(res, kind=c_int)
    call xcfun_eval_vec_C(fun, nr_points, density, d_pitch-1, &
         res, r_pitch-1)
  end subroutine
end module
