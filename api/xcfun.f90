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

  ! Evaluation modes
  enum, bind(C)
    enumerator :: XC_MODE_UNSET = 0
    enumerator XC_PARTIAL_DERIVATIVES
    enumerator XC_POTENTIAL
    enumerator XC_CONTRACTED
    enumerator XC_NR_MODES
  end enum

  ! Variable combinations
  enum, bind(C)
    enumerator :: XC_VARS_UNSET = -1
    enumerator XC_A
    enumerator XC_N
    enumerator XC_A_B
    enumerator XC_N_S
    enumerator XC_A_GAA
    enumerator XC_N_GNN
    enumerator XC_A_B_GAA_GAB_GBB
    enumerator XC_N_S_GNN_GNS_GSS
    enumerator XC_A_GAA_LAPA
    enumerator XC_A_GAA_TAUA
    enumerator XC_N_GNN_LAPN
    enumerator XC_N_GNN_TAUN
    enumerator XC_A_B_GAA_GAB_GBB_LAPA_LAPB
    enumerator XC_A_B_GAA_GAB_GBB_TAUA_TAUB
    enumerator XC_N_S_GNN_GNS_GSS_LAPN_LAPS
    enumerator XC_N_S_GNN_GNS_GSS_TAUN_TAUS
    enumerator XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB
    enumerator XC_A_B_GAA_GAB_GBB_LAPA_LAPB_TAUA_TAUB_JPAA_JPBB
    enumerator XC_N_S_GNN_GNS_GSS_LAPN_LAPS_TAUN_TAUS
    enumerator XC_A_AX_AY_AZ
    enumerator XC_A_B_AX_AY_AZ_BX_BY_BZ
    enumerator XC_N_NX_NY_NZ
    enumerator XC_N_S_NX_NY_NZ_SX_SY_SZ
    enumerator XC_A_AX_AY_AZ_TAUA
    enumerator XC_A_B_AX_AY_AZ_BX_BY_BZ_TAUA_TAUB
    enumerator XC_N_NX_NY_NZ_TAUN
    enumerator XC_N_S_NX_NY_NZ_SX_SY_SZ_TAUN_TAUS
    enumerator XC_A_2ND_TAYLOR
    enumerator XC_A_B_2ND_TAYLOR
    enumerator XC_N_2ND_TAYLOR
    enumerator XC_N_S_2ND_TAYLOR
  end enum

  private :: fstring_to_carray
  private :: text_handler

  interface
    function xcfun_version_C() result(text) &
      bind(C, name="xcfun_version")
      import
      type(c_ptr) :: text
    end function

    function xcfun_splash_C() result(text) &
      bind(C, name="xcfun_splash")
      import
      type(c_ptr) :: text
    end function

    function xcfun_authors_C() result(text) &
         bind(C, name="xcfun_authors")
      import
      type(c_ptr) :: text
    end function

    function xcfun_test() result(nfail) &
      bind(C)
      import
      integer(c_int) :: nfail
    end function

    function xcfun_is_compatible_library() result(is_compatible) &
      bind(C)
      import
      logical(c_bool) :: is_compatible
    end function

    function xcfun_which_vars(func_type, dens_type, laplacian, kinetic, current, explicit_derivatives) result(vars) &
         bind(C)
      import
      integer(c_int), intent(in), value :: func_type
      integer(c_int), intent(in), value :: dens_type
      integer(c_int), intent(in), value :: laplacian
      integer(c_int), intent(in), value :: kinetic
      integer(c_int), intent(in), value :: current
      integer(c_int), intent(in), value :: explicit_derivatives
      integer(c_int) :: vars
    end function

    function xcfun_which_mode(mode_type) result(mode) &
         bind(C)
      import
      integer(c_int), intent(in), value :: mode_type
      integer(c_int) :: mode
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
      type(c_ptr), intent(in), value :: fun
      character(kind=c_char, len=1), intent(in) :: name(*)
      real(c_double), intent(inout) :: val(*)
      integer(c_int) :: err
    end function

    function xcfun_is_gga(fun) result(is_gga) &
      bind(C)
      import
      type(c_ptr), intent(in), value :: fun
      logical(c_bool) :: is_gga
    end function

    function xcfun_is_metagga(fun) result(is_metagga) &
      bind(C)
      import
      type(c_ptr), intent(in), value :: fun
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

    subroutine xcfun_eval_vec_C(fun, nr_points, density, d_pitch, res, r_pitch) &
         bind(C, name="xcfun_eval_vec")
      import
      type(c_ptr), intent(in), value :: fun
      integer(c_int), intent(in), value :: nr_points
      real(c_double), intent(in) :: density(*)
      integer(c_int), intent(in), value :: d_pitch
      real(c_double), intent(inout) :: res(*)
      integer(c_int), intent(in), value :: r_pitch
    end subroutine xcfun_eval_vec_C
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

  function xcfun_version() result(text)
    character(kind=c_char, len=50) :: text

    text = text_handler(xcfun_version_C(), 50)
  end function

  function xcfun_splash() result(text)
    character(kind=c_char, len=5000) :: text

    text = text_handler(xcfun_splash_C(), 5000)
  end function

  function xcfun_authors() result(text)
    character(kind=c_char, len=5000) :: text

    text = text_handler(xcfun_authors_C(), 5000)
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
