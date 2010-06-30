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
  do n=1,XC_NR_PARAMS
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
! We use ilen and olen as the leading dimension (pitch) values,
! although with a single point these are not important.
  call xc_eval(funid,order,npoints,dens, ilen, funout, olen)
! Print output
  print *,'Output:'
  do n=1,olen
     print *,n,funout(n)
  enddo
end program
