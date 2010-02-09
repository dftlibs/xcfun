program xc_example

! this will calculate some functional derivatives
! for example densities

! this example contains all f90 interface toutines
! that are needed to "talk to" the xc_fun library

  use xc_fun

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

end program
