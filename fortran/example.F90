program xc_example

! this will calculate some functional derivatives
! for example densities

! this example contains all f90 interface toutines
! that are needed to "talk to" the xcfun library

  use xcfun
  integer funid,n,stat,ilen,olen,order
  character*1000 descr
  double precision, allocatable :: dens(:),funout(:)

  call xcfun_splash(descr)
  print *, descr(1:len_trim(descr))

! Create a new functional, we need this for interacting
! with the library.
  funid = xc_new_functional()

! Discover all settings
  print *,'All settings (functionals and parameters):'
  n = 0
  call xc_setting_name(funid,n,descr)
  do while(len_trim(descr).gt.0)
     print *,descr(1:len_trim(descr))
     n = n + 1
     call xc_setting_name(funid,n,descr)
  enddo

! Set up BLYP
  print *,'Setting up LDAERF'
!radovan: this works
! call xc_set_setting(funid,'beckex',1.0d0,stat)
! call xc_set_setting(funid,'lypc',1.0d0,stat)
!now trying:
  call xc_set_setting(funid, 'ldaerfx',    1.0d0, stat)
  call xc_set_setting(funid, 'ldaerfc',    1.0d0, stat)
!  call xc_set_setting(funid, 'lypc',       1.0d0, stat)

! Print currently set parameters
! Discover all settings
  print *,'Now defined (weights):'
  n = 0
  call xc_setting_name(funid,n,descr)
  do while(len_trim(descr).gt.0)
     if (xc_is_set(funid,descr)) then
        print *,descr(1:len_trim(descr)),xc_get_setting(funid,descr)
     endif
     n = n + 1
     call xc_setting_name(funid,n,descr)
  enddo
! Set alpha/beta density variable mode
  call xc_set_mode(funid,XC_VARS_AB)
! Type of functional
  print *,'Functional type:',xc_get_type(funid),' (0 LDA, 1 GGA, 2 MGGA)'
! Let's get derivatives to third order
  order = 2
  print *,'Order:',order
! Ask for the length of the input
  ilen = xc_input_length(funid)
  olen = xc_output_length(funid,order)
  print *,'Length of input:',ilen
  print *,'Length of Output:',olen
! Allocate data for evaluation
  allocate(dens(ilen))
  allocate(funout(olen))
! Set some random input
  do n=1,ilen
     dens(n) = n
  enddo
! Compute the functional and its derivatives up to order
  call xc_eval(funid,funout,order,dens)
! Print output
  print *,'Output:'
  do n=1,olen
     print *,n,funout(n)
  enddo

#if 0
  real(8)              :: fun(XC_PARAMS_LEN)
  real(8), allocatable :: derv(:)
  integer              :: max_order
  real(8)              :: n_a, n_b

  max_order = 2

! unset functional
  fun = 0.0d0

! set functional to pbe = pbex + pbec
  fun(XC_PBE_EXCHANGE)    = 1.0d0
  fun(XC_PBE_CORRELATION) = 1.0d0

! use ABFGH variable set
  call xc_set_functional(XC_ABFGH, fun)

! there are bin(2+5, 2) = 21 derivatives (including 0th) up to order 2
! allocate derv array to hold them
  allocate(derv(xc_len(XC_ABFGH, 2)))

!                         n_a
!                         |      n_b
!                         |      |      \nabla n_a \cdot \nabla n_b
!                         |      |      |      \nabla n_b \cdot \nabla n_b
!                         |      |      |      |      \nabla n_a \cdot \nabla n_b
!                         |      |      |      |      |
  call xc_eval(derv, 2, (/1.7d0, 1.7d0, 1.7d0, 1.7d0, 1.7d0/))

  print *, 'calculated', size(derv), 'PBE derivatives, for example:'
  print *, 'XC energy density: e = ', derv(1)
  print *, '(d/dn_a)^2 e         = ', derv(xc_index(XC_ABFGH, (/2, 0, 0, 0, 0/)))
  print *, 'd/dn_a d/dn_b e      = ', derv(xc_index(XC_ABFGH, (/1, 1, 0, 0, 0/)))

  n_a = 0.11273116096929683d0
  n_b = n_a

! reset functional to Slater
  fun = 0.0d0
  fun(XC_SLATER_EXCHANGE) = 1.0d0
  call xc_set_functional(XC_ABFGH, fun)
  call xc_eval(derv, 2, (/n_a, n_b, 0.0d0, 0.0d0, 0.0d0/))

  print *, 'some Slater derivatives:'
  print *, '(d/dn_a)^2 e         = ', derv(xc_index(XC_ABFGH, (/2, 0, 0, 0, 0/)))

! reset functional to VWN5
  fun = 0.0d0
  fun(XC_VWN5_CORRELATION) = 1.0d0
  call xc_set_functional(XC_ABFGH, fun)
  call xc_eval(derv, 2, (/n_a, n_b, 0.0d0, 0.0d0, 0.0d0/))

  print *, 'some VWN5 derivatives:'
  print *, '(d/dn_a)^2 e         = ', derv(xc_index(XC_ABFGH, (/2, 0, 0, 0, 0/)))

  deallocate(derv)
#endif
end program
