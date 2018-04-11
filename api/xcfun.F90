!
! XCFun, an arbitrary order exchange-correlation library
! Copyright (C) 2018 Ulf Ekström and contributors.
!
! This file is part of XCFun.
!
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
! For information on the complete list of contributors to the
! XCFun library, see: <https://xcfun.readthedocs.io/>
!

module xcfun
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

contains
  ! We pass strings as null terminated integer arrays to C, this
  ! should be portable if C and Fortran uses the same character set.
  subroutine str2ints(str, ints)
! ints must have length at least len(str)+1
    character, intent(in) :: str*(*)
    integer, intent(out) :: ints(*)
    integer i
    do i = 1, len_trim(str)
      ints(i) = ichar(str(i:i))
    enddo
    ints(len(str) + 1) = 0
  end subroutine str2ints

  subroutine ints2str(ints, str)
    character, intent(out) :: str*(*)
    integer, intent(in) :: ints(*)
    integer i, j
    i = 1
    do while (ints(i) .ne. 0)
      str(i:i) = char(ints(i))
      i = i + 1
    enddo
    do j = i, len(str)
      str(j:j) = ' '
    enddo
  end subroutine ints2str

  function xcfun_version()
    real(8) xcfun_version, xcfuve
    xcfun_version = xcfuve()
  end function xcfun_version

  subroutine xcfun_splash(text)
    character, intent(out) :: text*(*)
    integer :: ibuf(len(text) + 1)
    integer le
    le = len(text)
    call xcspla(ibuf, le)
    call ints2str(ibuf, text)
  end subroutine xcfun_splash

  function xcfun_sizeof_int()
    integer :: i, xcfun_sizeof_int
    if (huge(i) .eq. 2147483647) then
      xcfun_sizeof_int = 4
    else
      xcfun_sizeof_int = 8
    endif
  end function xcfun_sizeof_int

  ! Create a new, "empty" functional and return its id
  function xc_new_functional()
    integer :: xc_new_functional, xcnewf
    xc_new_functional = xcnewf(XCFUN_API_VERSION)
  end function xc_new_functional

  subroutine xc_free_functional(funid)
    integer funid
    call xcfree(funid)
  end subroutine xc_free_functional

  function xc_set(funid, param, val)
    character, intent(in) :: param*(*)
    integer, intent(in) :: funid
    integer :: xc_set, xcsets
    real(8), intent(in) :: val
    xc_set = xcsets(funid, val, len_trim(param), param)
  end function xc_set

  function xc_get(funid, param, val)
    character, intent(in) :: param*(*)
    integer, intent(in) :: funid
    integer :: xcgets, xc_get
    real(8), intent(out) :: val
    xc_get = xcgets(funid, val, len_trim(param), param)
  end function xc_get

  function xc_serialize(funid, buflen, buf)
    integer, intent(in) :: funid, buflen
    real(8), intent(out) :: buf(buflen)
    integer :: xc_serialize, xcseri
    xc_serialize = xcseri(funid, buflen, buf)
  end function xc_serialize

#ifdef XXX
  subroutine xc_short_description(param, description)
    integer, intent(in) :: param
    character, intent(out) :: description*(*)
    integer :: idescr(len(description) + 1)
    integer :: le
    le = len(description) + 1
    call xcssho(idescr, le, param)
    call ints2str(idescr, description)
  end subroutine

  subroutine xc_long_description(param, description)
    integer, intent(in) :: param
    character, intent(out) :: description*(*)
    integer :: idescr(len(description) + 1)
    integer :: le
    le = len(description) + 1
    call xcslon(idescr, le, param)
    call ints2str(idescr, description)
  end subroutine xc_long_description
#endif
  function xc_eval_setup(funid, vars, mode, order)
    integer :: funid, vars, mode, order, xc_eval_setup, xcevse
    xc_eval_setup = xcevse(funid, vars, mode, order)
  end function xc_eval_setup

  subroutine xc_eval(funid, npoints, densities, results)
    integer, intent(in) :: funid, npoints
    real(8), intent(in) :: densities(:, :)
    real(8), intent(out) :: results(:, :)

    integer :: i2

    i2 = 2
    IF (npoints .EQ. 1) i2 = 1

    call xceval(funid, npoints, densities(1, 1), densities(1, i2), &
                results(1, 1), results(1, i2))
  end subroutine xc_eval

  subroutine xc_eval_star(funid, npoints, m, densities, results)
    integer, intent(in) :: funid, npoints, m
    real(8), intent(in) :: densities(*)
    real(8), intent(out) :: results(:, :)
    call xceval(funid, npoints, densities(1), densities(m + 1), &
                results(1, 1), results(1, 2))
  end subroutine xc_eval_star

#if 0
  function xc_index(funid, exponents)
    integer, intent(in) :: exponents(*)
    integer funid, xc_index, xcdind
!radovan: it's easier to start in fortran with 1
!         it would be possible to start with 0
!         but then we would allocate 0:length-1
!         instead of 1:length
!         i find 0:length-1 allocations error-prone
!   xc_index = xcdind(funid,exponents)
    xc_index = xcdind(funid, exponents) + 1
  end function xc_index
#endif

  function xc_output_length(funid, order)
    integer, intent(in) :: funid, order
    integer xc_output_length, xcoule
    xc_output_length = xcoule(funid, order)
  end function xc_output_length
#if 0
  subroutine xc_param_name(setting_nr, name)
    integer, intent(in) :: setting_nr
    character, intent(out) :: name*(*)
    integer :: ibuf(len(name) + 1)
    integer :: le
    le = len(name) + 1
    call xcsnam(ibuf, le, setting_nr)
    call ints2str(ibuf, name)
  end subroutine xc_param_name
#endif
end module
