program xc_example

! this will calculate some functional derivatives
! for example densities

! this example contains all f90 interface toutines
! that are needed to "talk to" the xcfun library

  use xcfun
  integer funid,n,ilen,olen,order,npoints
  character*1000 descr
  double precision, allocatable :: dens(:),funout(:)

! Print some info and copyright about the library.
  call xcfun_splash(descr)
  print *, descr(1:len_trim(descr))

! Create a new functional, we need this for interacting
! with the library.
  funid = xc_new_functional()

! Set up BLYP
  print *,'Setting up BLYP'
  call xc_set_param(funid, XC_LYPC, 1.0d0)
  call xc_set_param(funid, XC_BECKEX, 1.0d0)

! Print currently set parameters
! Discover all settings
  print *,'Now defined (nonzero weights and parameters, including defaults):'
  do n=0,XC_NR_PARAMS-1
     call xc_param_name(n,descr)
     if (xc_get_param(funid,n).ne.0.0D0) then
        print *, descr(4:len_trim(descr)),xc_get_param(funid,n)
     endif
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
  npoints = 1 ! We have only one point in this example
! Compute the functional and its derivatives up to order
  call xc_eval(funid,order,npoints,dens, funout)
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
