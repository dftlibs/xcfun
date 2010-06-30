
module xcfun
  use xcfun_autogen
  implicit none
! These parameters should mirror those in xcfun.h
  integer, parameter ::  XC_VARS_A  = 0
  integer, parameter ::  XC_VARS_R  = 1
  integer, parameter ::  XC_VARS_AB = 2
  integer, parameter ::  XC_VARS_RS = 3
  integer, parameter ::  XC_LDA  = 0
  integer, parameter ::  XC_GGA  = 1
  integer, parameter ::  XC_MGGA = 2
! Indices into the output array of derivatives. Fortran numbering.
  integer, parameter :: XC_D0 = 1
  integer, parameter :: XC_D1 = 2
  integer, parameter :: XC_D2 = 3

! lda derivatives, all odd spin density orders give zero at closed-shell reference
  integer, parameter :: XC_D00 =  1
  integer, parameter :: XC_D10 =  2
  integer, parameter :: XC_D01 =  3 !zero at closed-shell reference with XC_VARS_RS 
  integer, parameter :: XC_D20 =  4
  integer, parameter :: XC_D11 =  5 !zero at closed-shell reference with XC_VARS_RS 
  integer, parameter :: XC_D02 =  6
  integer, parameter :: XC_D30 =  7
  integer, parameter :: XC_D21 =  8 !zero at closed-shell reference with XC_VARS_RS 
  integer, parameter :: XC_D12 =  9
  integer, parameter :: XC_D03 = 10 !zero at closed-shell reference with XC_VARS_RS 

  integer, parameter :: XC_D00000 = 1
  integer, parameter :: XC_D10000 = 2
  integer, parameter :: XC_D01000 = 3
  integer, parameter :: XC_D00100 = 4
  integer, parameter :: XC_D00010 = 5
  integer, parameter :: XC_D00001 = 6
  integer, parameter :: XC_D20000 = 7
  integer, parameter :: XC_D11000 = 8
  integer, parameter :: XC_D10100 = 9
  integer, parameter :: XC_D10010 = 10
  integer, parameter :: XC_D10001 = 11
  integer, parameter :: XC_D02000 = 12
  integer, parameter :: XC_D01100 = 13
  integer, parameter :: XC_D01010 = 14
  integer, parameter :: XC_D01001 = 15
  integer, parameter :: XC_D00200 = 16
  integer, parameter :: XC_D00110 = 17
  integer, parameter :: XC_D00101 = 18
  integer, parameter :: XC_D00020 = 19
  integer, parameter :: XC_D00011 = 20
  integer, parameter :: XC_D00002 = 21

!radovan: here are the pointers to 3rd gga derivatives
!         be careful, only those that are nonzero at closed-shell reference 
!         with XC_VARS_RS are here
  integer, parameter :: XC_D30000 = 22
  integer, parameter :: XC_D20100 = 24
  integer, parameter :: XC_D12000 = 27
  integer, parameter :: XC_D11001 = 30
  integer, parameter :: XC_D10200 = 31
  integer, parameter :: XC_D10002 = 36
  integer, parameter :: XC_D02100 = 38
  integer, parameter :: XC_D01101 = 43
  integer, parameter :: XC_D00300 = 47
  integer, parameter :: XC_D00102 = 52

contains
  ! We pass strings as null terminated integer arrays to C, this
  ! should be portable if C and Fortran uses the same character set.
  subroutine str2ints(str,ints)
! ints must have length at least len(str)+1
    character, intent(in) :: str*(*)
    integer, intent(out) :: ints(*)
    integer i
    do i=1,len_trim(str)
       ints(i) = ichar(str(i:i))
    enddo
    ints(len(str)+1) = 0
  end subroutine str2ints

  subroutine ints2str(ints,str)
    character, intent(out) :: str*(*)
    integer, intent(in) :: ints(*)
    integer i,j
    i = 1
    do while (ints(i).ne.0)
       str(i:i) = char(ints(i))
       i = i + 1
    enddo
    do j=i,len(str)
       str(j:j) = ' '
    enddo
  end subroutine ints2str


  function xcfun_version()
    double precision xcfun_version, xcfuve
    xcfun_version = xcfuve()
  end function xcfun_version

  subroutine xcfun_splash(text)
    character, intent(out) :: text*(*)
    integer :: ibuf(len(text)+1)
    integer le
    le = len(text)
    call xcspla(ibuf,le)
    call ints2str(ibuf,text)
  end subroutine xcfun_splash

  ! Create a new, "empty" functional and return its id
  function xc_new_functional()
    integer xc_new_functional, xcnewf
    xc_new_functional = xcnewf()
  end function xc_new_functional

  subroutine xc_free_functional(funid)
    integer funid
    call xcfree(funid)
  end subroutine xc_free_functional

  subroutine xc_set_param(funid, param, val)
    integer, intent(in) :: funid, param
    double precision, intent(in) :: val
    call xcsets(funid,param,val)
  end subroutine xc_set_param
  
  function xc_get_param(funid, param)
    integer, intent(in) :: funid, param
    double precision xcgets, xc_get_param
    xc_get_param = xcgets(funid,param)
  end function xc_get_param

  function xc_is_functional(param)
    logical :: xc_is_functional
    integer, intent(in) :: param
    integer xcisfu,res
    res = xcisfu(param)
    if (res.eq.0) then
       xc_is_functional = .false.
    else
       xc_is_functional = .true.
    endif
  end function

  subroutine xc_short_description(param,description)
    integer, intent(in) :: param
    character, intent(out) :: description*(*)
    integer :: idescr(len(description)+1)
    integer :: le
    le = len(description)+1
    call xcssho(idescr,le,param)
    call ints2str(idescr,description)
  end subroutine

  subroutine xc_long_description(param,description)
    integer, intent(in) :: param
    character, intent(out) :: description*(*)
    integer :: idescr(len(description)+1)
    integer :: le
    le = len(description)+1
    call xcslon(idescr,le,param)
    call ints2str(idescr,description)
  end subroutine

  subroutine xc_set_mode(funid, mode)
    integer, intent(in) :: funid, mode
    call xcsmod(funid,mode)
  end subroutine xc_set_mode

  function xc_get_type(funid)
    integer, intent(in) :: funid
    integer xc_get_type, xcgett
    xc_get_type = xcgett(funid)
  end function xc_get_type

  subroutine xc_eval(funid, order, npoints, densvars, denspitch, res, respitch)
    integer, intent(in) :: funid, npoints, order,denspitch,respitch
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
    call xceval(funid,order,npoints,densvars,denspitch,res,respitch)
  end subroutine xc_eval

  function xc_index(funid, exponents)
    integer, intent(in) :: exponents(*)
    integer funid,xc_index, xcdind
!radovan: it's easier to start in fortran with 1
!         it would be possible to start with 0
!         but then we would allocate 0:length-1
!         instead of 1:length
!         i find 0:length-1 allocations error-prone
!   xc_index = xcdind(funid,exponents)
    xc_index = xcdind(funid,exponents) + 1
  end function xc_index

  function xc_input_length(funid)
    integer, intent(in) :: funid
    integer xc_input_length, xcinle
    xc_input_length = xcinle(funid)
  end function xc_input_length

  function xc_output_length(funid, order)
    integer, intent(in) :: funid, order
    integer xc_output_length, xcoule
    xc_output_length = xcoule(funid,order)
  end function xc_output_length

  subroutine xc_param_name(setting_nr, name)
    integer, intent(in) :: setting_nr
    character, intent(out) :: name*(*)
    integer :: ibuf(len(name)+1)
    integer :: le
    le = len(name)+1
    call xcsnam(ibuf,le,setting_nr)
    call ints2str(ibuf,name)
  end subroutine xc_param_name
end module
