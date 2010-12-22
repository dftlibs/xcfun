
module xcfun
  use xcfun_autogen
  implicit none

! these parameters should mirror those in xcfun.h
  integer, parameter ::  XC_VARS_A  = 0
  integer, parameter ::  XC_VARS_N  = 1
  integer, parameter ::  XC_VARS_AB = 2
  integer, parameter ::  XC_VARS_NS = 3
  integer, parameter ::  XC_LDA  = 0
  integer, parameter ::  XC_GGA  = 1
  integer, parameter ::  XC_MGGA = 2

! one-variable derivatives up to order 4
  integer, parameter :: XC_D0 = 1
  integer, parameter :: XC_D1 = 2
  integer, parameter :: XC_D2 = 3
  integer, parameter :: XC_D3 = 4
  integer, parameter :: XC_D4 = 5

! two-variable derivatives up to order 4
  integer, parameter :: XC_D00 =  1
  integer, parameter :: XC_D10 =  2
  integer, parameter :: XC_D01 =  3
  integer, parameter :: XC_D20 =  4
  integer, parameter :: XC_D11 =  5
  integer, parameter :: XC_D02 =  6
  integer, parameter :: XC_D30 =  7
  integer, parameter :: XC_D21 =  8
  integer, parameter :: XC_D12 =  9
  integer, parameter :: XC_D03 = 10
  integer, parameter :: XC_D40 = 11
  integer, parameter :: XC_D31 = 12
  integer, parameter :: XC_D22 = 13
  integer, parameter :: XC_D13 = 14
  integer, parameter :: XC_D04 = 15

! five-variable derivatives up to order 4
! all up to order 3 are here
! not all order 4
  integer, parameter :: XC_D00000 =   1
  integer, parameter :: XC_D10000 =   2
  integer, parameter :: XC_D01000 =   3
  integer, parameter :: XC_D00100 =   4
  integer, parameter :: XC_D00010 =   5
  integer, parameter :: XC_D00001 =   6
  integer, parameter :: XC_D20000 =   7
  integer, parameter :: XC_D11000 =   8
  integer, parameter :: XC_D10100 =   9
  integer, parameter :: XC_D10010 =  10
  integer, parameter :: XC_D10001 =  11
  integer, parameter :: XC_D02000 =  12
  integer, parameter :: XC_D01100 =  13
  integer, parameter :: XC_D01010 =  14
  integer, parameter :: XC_D01001 =  15
  integer, parameter :: XC_D00200 =  16
  integer, parameter :: XC_D00110 =  17
  integer, parameter :: XC_D00101 =  18
  integer, parameter :: XC_D00020 =  19
  integer, parameter :: XC_D00011 =  20
  integer, parameter :: XC_D00002 =  21
  integer, parameter :: XC_D30000 =  22
  integer, parameter :: XC_D21000 =  23
  integer, parameter :: XC_D20100 =  24
  integer, parameter :: XC_D20010 =  25
  integer, parameter :: XC_D20001 =  26
  integer, parameter :: XC_D12000 =  27
  integer, parameter :: XC_D11100 =  28
  integer, parameter :: XC_D11010 =  29
  integer, parameter :: XC_D11001 =  30
  integer, parameter :: XC_D10200 =  31
  integer, parameter :: XC_D10110 =  32
  integer, parameter :: XC_D10101 =  33
  integer, parameter :: XC_D10020 =  34
  integer, parameter :: XC_D10011 =  35
  integer, parameter :: XC_D10002 =  36
  integer, parameter :: XC_D03000 =  37
  integer, parameter :: XC_D02100 =  38
  integer, parameter :: XC_D02010 =  39
  integer, parameter :: XC_D02001 =  40
  integer, parameter :: XC_D01200 =  41
  integer, parameter :: XC_D01110 =  42
  integer, parameter :: XC_D01101 =  43
  integer, parameter :: XC_D01020 =  44
  integer, parameter :: XC_D01011 =  45
  integer, parameter :: XC_D01002 =  46
  integer, parameter :: XC_D00300 =  47
  integer, parameter :: XC_D00210 =  48
  integer, parameter :: XC_D00201 =  49
  integer, parameter :: XC_D00120 =  50
  integer, parameter :: XC_D00111 =  51
  integer, parameter :: XC_D00102 =  52
  integer, parameter :: XC_D00030 =  53
  integer, parameter :: XC_D00021 =  54
  integer, parameter :: XC_D00012 =  55
  integer, parameter :: XC_D00003 =  56
  integer, parameter :: XC_D40000 =  57
  integer, parameter :: XC_D30100 =  59
  integer, parameter :: XC_D20200 =  66
  integer, parameter :: XC_D12100 =  73
  integer, parameter :: XC_D12010 =  74
  integer, parameter :: XC_D12001 =  75
  integer, parameter :: XC_D10300 =  82
  integer, parameter :: XC_D00400 = 112

! seven-variable derivatives up to order 4
! so far only linear response
  integer, parameter :: XC_D0000000 =  1
  integer, parameter :: XC_D1000000 =  2
  integer, parameter :: XC_D0010000 =  4
  integer, parameter :: XC_D0000100 =  6
  integer, parameter :: XC_D0000010 =  7
  integer, parameter :: XC_D0000001 =  8
  integer, parameter :: XC_D2000000 =  9
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

  subroutine xc_eval(funid, order, npoints, densities, results)
    integer, intent(in) :: funid, order, npoints
    double precision, intent(in) :: densities(:,:)
    double precision, intent(out) :: results(:,:)
    call xceval(funid,order,npoints,densities(1,1),densities(1,2),&
         results(1,1),results(1,2))
  end subroutine xc_eval

  subroutine xc_potential(funid, density, energy, potentials)
    integer, intent(in) :: funid    
    double precision, intent(in) :: density(:)
    double precision, intent(out) :: energy, potentials(:)
    call xcpotential(funid,density,energy,potentials)
  end subroutine xc_potential

  function xc_index(funid, exponents)
    integer, intent(in) :: exponents(*)
    integer funid,xc_index, xcdind
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
