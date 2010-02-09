!The following variables are used in the code:
! A = \rho_\alpha
! B = \rho_\beta
! F = \nabla A \cdot \nabla A (usually called \sigma_\alpha\alpha)
! G = \nabla B \cdot \nabla B (\sigma_\beta\beta)
! H = \nabla A \cdot \nabla B (\sigma_\alpha\beta)
!Alternatively you can use;
! R = \rho_\alpha + \rho_\beta
! S = \rho_\alpha - \rho_\beta
! Z = \nabla R \cdot \nabla R (sometimes called \sigma)
! X = \nabla S \cdot \nabla S
! Y = \nabla R \cdot \nabla S

module xc_fun
  implicit none
  integer, parameter ::  XC_MAX_ORDER = 6

! Index of each functional parameter in the argument
! of set_functional
  integer, parameter ::  XC_VWN5_CORRELATION = 1
  integer, parameter ::  XC_SLATER_EXCHANGE  = 2
  integer, parameter ::  XC_LYP_CORRELATION  = 3
  integer, parameter ::  XC_PW91_EXCHANGE    = 4
  integer, parameter ::  XC_PBE_CORRELATION  = 5
  integer, parameter ::  XC_OLD_PBE_CORRELATION = 6
  integer, parameter ::  XC_BECKE_EXCHANGE   = 7
  integer, parameter ::  XC_OLD_PBE_EXCHANGE = 8
  integer, parameter ::  XC_OLD_RPBE_EXCHANGE= 9 
  integer, parameter ::  XC_PBE_EXCHANGE     = 10
  integer, parameter ::  XC_REVPBE_EXCHANGE  = 11
  integer, parameter ::  XC_PW92_CORRELATION = 12
  integer, parameter ::  KIN_TF              = 13
  integer, parameter ::  KIN_PW91            = 14
! Make sure XC_PARAMS_LEN matches the number of parameters above
  integer, parameter ::  XC_PARAMS_LEN = 14 

! The different modes (density variables) xc_fun can work with,
! these are passed as the first argument to set_functional()
   integer, parameter ::   XC_A  = 0 !Fully spin polarized LDA
   integer, parameter ::   XC_R  = 1 !Unpolarized LDA
   integer, parameter ::   XC_AB = 2 !LDA with alpha/beta densities
   integer, parameter ::   XC_RS = 3 !LDA with total density/polarization
   integer, parameter ::   XC_AF = 4 !Fully spin polarized GGA
   integer, parameter ::   XC_RZ = 5 !Unpolarized GGA
   integer, parameter ::   XC_ABFGH = 6 !Generic GGA alpha,beta variables 
   integer, parameter ::   XC_RSZXY = 7 !Generic GGA R,S variables
   integer, parameter ::   XC_R_ZI_ZJ_ZK = 8 !Generic GGA R,S variables	

contains
! Set up the functional used by xc_eval. The mode is one of 
! XC_A,XC_R, .. defined above. It defines the variables that
! will be used as input and output to xc_eval. For example
! using XC_ABFGH will require densvars to be an array of
! the five variables A,B,F,G,H described above. The output
! derivatives are given wrt these same variables. 
  subroutine xc_set_functional(mode,params)
    integer, intent(in) :: mode
    double precision, intent(in) :: params(XC_PARAMS_LEN)
    call xc_set_functional_fortran(mode,params,XC_PARAMS_LEN)
  end subroutine xc_set_functional

! Eval the functional previously set by xc_set_functional.
! The size and layout of res and densvars depends on the
! mode given to xc_set_functional.
  subroutine xc_eval(res, order, densvars)
    integer order
    double precision, intent(out) :: res(*)
    double precision, intent(in) :: densvars(*)
    call xc_eval_fortran(res,order,densvars)
  end subroutine xc_eval

  function xc_index(mode, exponents)
    integer, intent(in) :: mode
    integer, intent(in) :: exponents(*)
    integer xc_index, xc_index_fortran
    xc_index = xc_index_fortran(mode,exponents)
  end function xc_index

  function xc_len(mode, order)
    integer, intent(in) :: mode
    integer, intent(in) :: order
    integer xc_len, xc_len_fortran
    xc_len = xc_len_fortran(mode,order)
  end function xc_len
end module
