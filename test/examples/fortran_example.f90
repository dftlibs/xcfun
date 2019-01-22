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


subroutine assert(predicate, error_message)

  implicit none

  logical, intent(in) :: predicate
  character(*), intent(in) :: error_message

  if (.not. predicate) then
    print *, 'ERROR:', error_message
    stop 1
  endif

end subroutine


program xc_example

! This example contains calls to f90 interface routines that are needed to "talk
! to" the xcfun library and demonstrates how to use them.

! We will compute the kernel for an unpolarized system using total density and
! the gradient components as the variables. These are linear in the density
! matrix, which helps the code using the results from xcfun.

  use xcfun, only: XC_CONTRACTED, &
                   XC_N_NX_NY_NZ, &
                   xcfun_splash, &
                   xc_new_functional, &
                   xc_set, &
                   xc_eval_setup, &
                   xc_eval, &
                   xc_free_functional

  implicit none

  integer, parameter :: num_points = 1
  integer, parameter :: num_variables = 4  ! we will use XC_N_NX_NY_NZ

  character(1000) :: text
  integer :: id, order, ierr, ideriv, ipoint
  integer :: vector_length
  real(8) :: res

  real(8), allocatable :: density(:, :, :)
  real(8), allocatable :: input_array(:, :)
  real(8), allocatable :: output_array(:, :)
  real(8), allocatable :: res_reference(:)


  ! print some info and copyright about the library
  ! please always include this info in your code
  call xcfun_splash(text)
  print *, text(1:len_trim(text))

  ! create a new functional
  ! we need this for interacting with the library
  id = xc_new_functional()

  ! in this example we use PBE
  print *, 'Setting up PBE'
  ierr = xc_set(id, 'pbe', 1.0d0)
  call assert(ierr == 0, "functional name not recognized")

  ! First we just compute the energy, i.e. the 0-th order integrand.
  ! We have one gridpoint, and four variables, density N and gradient
  ! components NX NY NZ.
  order = 0
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate (density(vector_length, num_variables, num_points))
  density = 0.0d0
  do ipoint = 1, num_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 !        n
    density(1, 2, ipoint) = 2.0d0 !nabla_x n
    density(1, 3, ipoint) = 3.0d0 !nabla_y n
    density(1, 4, ipoint) = 4.0d0 !nabla_z n
  end do

  allocate (input_array(num_variables*vector_length, num_points))
  allocate (output_array(vector_length, num_points))

  do ipoint = 1, num_points
    input_array(:, ipoint) = (/density(:, 1, ipoint), &
                               density(:, 2, ipoint), &
                               density(:, 3, ipoint), &
                               density(:, 4, ipoint)/)
  end do
  call xc_eval(id, num_points, input_array, output_array)
  res = output_array(vector_length, 1)
  print *, 'The XC energy density is', res

  deallocate (density)
  deallocate (input_array)
  deallocate (output_array)

  ! now let's compute the first derivatives ('potential')
  ! and contract them with the first order densities
  ! first set up xcfun for first derivatives
  order = 1
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate (density(vector_length, num_variables, num_points))
  density = 0.0d0
  do ipoint = 1, num_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 !        n zeroth order
    density(1, 2, ipoint) = 2.0d0 !nabla_x n zeroth order
    density(1, 3, ipoint) = 3.0d0 !nabla_y n zeroth order
    density(1, 4, ipoint) = 4.0d0 !nabla_z n zeroth order
    density(2, 1, ipoint) = 5.0d0 !        n first order
    density(2, 2, ipoint) = 6.0d0 !nabla_x n first order
    density(2, 3, ipoint) = 7.0d0 !nabla_y n first order
    density(2, 4, ipoint) = 8.0d0 !nabla_z n first order
  end do

  allocate (input_array(num_variables*vector_length, num_points))
  allocate (output_array(vector_length, num_points))

!     It is possible to put in other numbers than 1 and 0
!     this is where perturbed densities go for automatic contraction

    do ipoint = 1, num_points
      input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                 density(:, 2, ipoint), &
                                 density(:, 3, ipoint), &
                                 density(:, 4, ipoint)/)
    end do
    call xc_eval(id, num_points, input_array, output_array)

  deallocate (density)
  deallocate (input_array)
  deallocate (output_array)


  ! reference values for self test
  allocate (res_reference(num_variables))
  res_reference(1) = -0.97182634532897016d0
  res_reference(2) = -8.9803044914016916d-3
  res_reference(3) = -1.3470456737102541d-2
  res_reference(4) = -1.7960608982803383d-2


  ! now let's compute the first derivatives ('potential')
  ! and do NOT contract them!
  ! first set up xcfun for first derivatives
  order = 1
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate (density(vector_length, num_variables, num_points))
  density = 0.0d0
  do ipoint = 1, num_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 !        n zeroth order
    density(1, 2, ipoint) = 2.0d0 !nabla_x n zeroth order
    density(1, 3, ipoint) = 3.0d0 !nabla_y n zeroth order
    density(1, 4, ipoint) = 4.0d0 !nabla_z n zeroth order
  end do

  allocate (input_array(num_variables*vector_length, num_points))
  allocate (output_array(vector_length, num_points))

  ! here is a tricky part, we need a loop and we need to put
  ! in a '1' where we want the derivative
  ! we loop over all variables we want a derivative with respect to
  do ideriv = 1, num_variables
    do ipoint = 1, num_points
      density(2, :, ipoint) = 0.0d0
      density(2, ideriv, ipoint) = 1.0d0
    end do
    ! It is possible to put in other numbers than 1 and 0
    ! this is where perturbed densities go for automatic contraction

    do ipoint = 1, num_points
      input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                 density(:, 2, ipoint), &
                                 density(:, 3, ipoint), &
                                 density(:, 4, ipoint)/)
    end do
    call xc_eval(id, num_points, input_array, output_array)

    res = output_array(vector_length, 1)
    print *, 'dE/dx_i for i =', ideriv, 'is', res

    ! test against reference numbers
    call assert((abs(res - res_reference(ideriv)) < 1.0d-6), &
                "derivatives do not match reference numbers")
  end do

  deallocate (density)
  deallocate (input_array)
  deallocate (output_array)
  deallocate (res_reference)

  ! reference values for self test
  allocate (res_reference(num_variables))
  res_reference(1) = -2.0795476275887927d0
  res_reference(2) = 2.1489169824673575d-2
  res_reference(3) = 4.1214059228412030d-2
  res_reference(4) = 6.0938948632150471d-2


  ! now second derivative of the "potential", contracted with ONE perturbed density.
  ! hopefully the strange input/output_array format will start to make sense
  order = 2
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate (density(vector_length, num_variables, num_points))
  density = 0.0d0
  do ipoint = 1, num_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 ! zeroth order
    density(1, 2, ipoint) = 2.0d0 ! zeroth order
    density(1, 3, ipoint) = 3.0d0 ! zeroth order
    density(1, 4, ipoint) = 4.0d0 ! zeroth order
    density(2, 1, ipoint) = 5.0d0 ! first order
    density(2, 2, ipoint) = 6.0d0 ! first order
    density(2, 3, ipoint) = 7.0d0 ! first order
    density(2, 4, ipoint) = 8.0d0 ! first order
  end do

  allocate (input_array(num_variables*vector_length, num_points))
  allocate (output_array(vector_length, num_points))

  ! now we have to put in ones and zeroes again to generate derivatives (which then generate matrix elements in ADF)
  ! density(1, :, :)  unperturbed values
  ! density(2, :, :)  perturbed values
  ! density(3, :, :)  ones and zeros
  ! density(4, :, :)  not used in this case
  ! For the second order case, there will be 2^order elements, 1 ground state, 2^order-2 perturbed and one unused.
  do ideriv = 1, num_variables
    do ipoint = 1, num_points
      density(3, :, ipoint) = 0.0d0
      density(3, ideriv, ipoint) = 1.0d0
    end do
    ! It is possible to put in other numbers than 1 and 0 in the density2 array
    ! this is where perturbed densities go for automatic contraction

    do ipoint = 1, num_points
      input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                 density(:, 2, ipoint), &
                                 density(:, 3, ipoint), &
                                 density(:, 4, ipoint)/)
    end do
    call xc_eval(id, num_points, input_array, output_array)

    res = output_array(vector_length, 1)
    print *, 'd^2 E/dx_i dD_pert for i =', ideriv, 'is', res

    ! test against reference numbers
    call assert((abs(res - res_reference(ideriv)) < 1.0d-6), &
                "derivatives do not match reference numbers")
  end do

  deallocate (density)
  deallocate (input_array)
  deallocate (output_array)
  deallocate (res_reference)

!  The computed second derivatives are now already contracted with the
!  perturbed density! This works to any order.  In this example we put in ones and
!  zeros, which you need to do to compute matrix elements. You can also put in a
!  perturbed density instead of ones and zeros, then you generate a response
!  contribution d^2E/dD1dD2. Typically this is what you want to do to high order,
!  because by the 2N+1 rule you don't have to compute matrix elements to very high
!  order. Then you need only one call of xc_eval, not one for each variable. This
!  is very efficient, much better than computing all partial derivatives. For
!  second derivatives it might be better to compute all partial derivatives and
!  reuse them in the response solver.
!
!  In this example we use the variables that are linear in the density matrix.
!  This makes life easier because you can construct them trivially from perturbed
!  density matrices. If you use i.e. the square norm of the density gradient
!  things get more complicated.
!
!  For the spin polarized cases (using density and spin density or up- and down
!  densities) we would use 8 variables instead of 4.
!
!  This will be exemplified in the following.



  ! Now contract the 2nd derivative completely
  !  allocate density and fill in some fantasy values
  !  Now the shape of zdensity is different from the last example.
  !  Typically, we want to contract the first kernel with e.g.
  !  first and/or second-order perturbed densities and again, we have
  !  2^order, i.e. 4 elements:
  !  1.) ground state density (zeroth order)
  !  2.) perturbed density a) (first order)
  !  3.) perturbed density b) (first order)
  !  4.) perturbed density ab) (second order)
  allocate (density(vector_length, num_variables, num_points))
  density = 0.0d0
  do ipoint = 1, num_points
    density(1, 1, ipoint) = 1.0d0    ! zeroth order
    density(1, 2, ipoint) = 2.0d0    ! zeroth order
    density(1, 3, ipoint) = 3.0d0    ! zeroth order
    density(1, 4, ipoint) = 4.0d0    ! zeroth order
    density(2, 1, ipoint) = 5.0d0    ! first order
    density(2, 2, ipoint) = 6.0d0   ! first order
    density(2, 3, ipoint) = 7.0d0   ! first order
    density(2, 4, ipoint) = 8.0d0   ! first order
    density(3, 1, ipoint) = 5.0d0    ! first order
    density(3, 2, ipoint) = 6.0d0    ! first order
    density(3, 3, ipoint) = 7.0d0    ! first order
    density(3, 4, ipoint) = 8.0d0    ! first order
    density(4, 1, ipoint) = 0.0d0    ! first order
    density(4, 2, ipoint) = 0.0d0    ! first order
    density(4, 3, ipoint) = 0.0d0    ! first order
    density(4, 4, ipoint) = 0.0d0    ! first order
  end do

  allocate (input_array(num_variables*vector_length, num_points))
  allocate (output_array(vector_length, num_points))

  do ipoint = 1, num_points
    input_array(:, ipoint) = (/density(:, 1, ipoint), &
                               density(:, 2, ipoint), &
                               density(:, 3, ipoint), &
                               density(:, 4, ipoint)/)
  end do
  call xc_eval(id, num_points, input_array, output_array)

  !  The output array will then contain:
  ! 1.) Product of functional with ground state density
  ! 2.) Product of first-order functional derivative with first-order density a)
  ! 3.) Product of first-order functional derivative with first-order density b)
  ! 4.) Product of first-order functional derivative with second-order density ab) +
  !     Product of second-order functional derivative with first-order densities a) and b)

  deallocate (density)
  deallocate (input_array)
  deallocate (output_array)

  !  now third derivative of the "potential", contracted with perturbed densities.
  order = 3
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  !  allocate density and fill in some fantasy values
  allocate (density(vector_length, num_variables, num_points))
  density = 0.0d0
  do ipoint = 1, num_points
    density(1, 1, ipoint) = 1.0d0   ! zeroth order
    density(1, 2, ipoint) = 2.0d0   ! zeroth order
    density(1, 3, ipoint) = 3.0d0   ! zeroth order
    density(1, 4, ipoint) = 4.0d0   ! zeroth order
    density(2, 1, ipoint) = 5.0d0   ! first order (1)
    density(2, 2, ipoint) = 6.0d0   ! first order (1)
    density(2, 3, ipoint) = 7.0d0   ! first order (1)
    density(2, 4, ipoint) = 8.0d0   ! first order (1)
    density(3, 1, ipoint) = 9.0d0   ! first order (2)
    density(3, 2, ipoint) = 10.0d0  ! first order (2)
    density(3, 3, ipoint) = 11.0d0  ! first order (2)
    density(3, 4, ipoint) = 12.0d0  ! first order (2)
    density(4, 1, ipoint) = 5.0d0   ! second order (depending on (1) and (2))
    density(4, 2, ipoint) = 6.0d0   ! second order (depending on (1) and (2))
    density(4, 3, ipoint) = 7.0d0   ! second order (depending on (1) and (2))
    density(4, 4, ipoint) = 8.0d0   ! second order (depending on (1) and (2))
    density(5, 1, ipoint) = 9.0d0   ! first order (3)
    density(5, 2, ipoint) = 10.0d0  ! first order (3)
    density(5, 3, ipoint) = 11.0d0  ! first order (3)
    density(5, 4, ipoint) = 12.0d0  ! first order (3)
    density(6, 1, ipoint) = 5.0d0   ! second order (depending on (1) and (3))
    density(6, 2, ipoint) = 6.0d0   ! second order (depending on (1) and (3))
    density(6, 3, ipoint) = 7.0d0   ! second order (depending on (1) and (3))
    density(6, 4, ipoint) = 8.0d0   ! second order (depending on (1) and (3))
    density(7, 1, ipoint) = 9.0d0   ! second order (depending on (2) and (3))
    density(7, 2, ipoint) = 10.0d0  ! second order (depending on (2) and (3))
    density(7, 3, ipoint) = 11.0d0  ! second order (depending on (2) and (3))
    density(7, 4, ipoint) = 12.0d0  ! second order (depending on (2) and (3))
    density(8, 1, ipoint) = 0.0d0   ! third order (depending on (1), (2) and (3))
    density(8, 2, ipoint) = 0.0d0   ! third order (depending on (1), (2) and (3))
    density(8, 3, ipoint) = 0.0d0   ! third order (depending on (1), (2) and (3))
    density(8, 4, ipoint) = 0.0d0   ! third order (depending on (1), (2) and (3))
  end do

  allocate (input_array(num_variables*vector_length, num_points))
  allocate (output_array(vector_length, num_points))

    do ipoint = 1, num_points
      input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                 density(:, 2, ipoint), &
                                 density(:, 3, ipoint), &
                                 density(:, 4, ipoint)/)
    end do

    call xc_eval(id, num_points, input_array, output_array)


  deallocate (density)
  deallocate (input_array)
  deallocate (output_array)

  !Release the functional
  call xc_free_functional(id)

  print *, 'Kernel test has ended properly!'

end program
