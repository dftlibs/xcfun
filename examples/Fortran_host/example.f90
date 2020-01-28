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

  ! we consider only one grid point
  integer, parameter :: num_grid_points = 1

  ! we will use XC_N_NX_NY_NZ
  ! N: density
  ! NX: x-gradient of the density
  ! NY: y-gradient of the density
  ! NZ: z-gradient of the density
  integer, parameter :: num_density_variables = 4

  character(1000) :: text
  integer :: id, order, ierr, ipoint
  integer :: vector_length
  real(8) :: res

  real(8), allocatable :: density(:, :, :)


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


  !-----------------------------------------------------------------------------
  ! in the first example we compute the XC energy ("order 0 derivative")

  order = 0
  ! XC_CONTRACTED here has nothing to do with contracted basis sets
  ! it means we will evaluated in the XC_CONTRACTED mode and
  ! internally contract functional derivatives with the density taylor expansion
  ! in other words: we will not have to explicitly assemble/contract partial
  ! derivatives outside of XCFun
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate(density(vector_length, num_density_variables, num_grid_points))
  density = 0.0d0
  do ipoint = 1, num_grid_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 !         n
    density(1, 2, ipoint) = 2.0d0 ! nabla_x n
    density(1, 3, ipoint) = 3.0d0 ! nabla_y n
    density(1, 4, ipoint) = 4.0d0 ! nabla_z n
  end do

  res = derivative(id, &
                   num_grid_points, &
                   num_density_variables, &
                   vector_length, &
                   density)
  deallocate(density)
  print *, 'The XC energy density is', res

  ! compare with reference
  call assert((abs(-0.86494159400066051d0 - res) < 1.0d-6), &
              "derivatives do not match reference numbers")


  !-----------------------------------------------------------------------------
  ! now we will compute the first derivatives ('potential')
  ! and contract them with the first order densities

  order = 1
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate(density(vector_length, num_density_variables, num_grid_points))
  density = 0.0d0
  do ipoint = 1, num_grid_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 !         n zeroth order
    density(1, 2, ipoint) = 2.0d0 ! nabla_x n zeroth order
    density(1, 3, ipoint) = 3.0d0 ! nabla_y n zeroth order
    density(1, 4, ipoint) = 4.0d0 ! nabla_z n zeroth order
    density(2, 1, ipoint) = 5.0d0 !         n first order
    density(2, 2, ipoint) = 6.0d0 ! nabla_x n first order
    density(2, 3, ipoint) = 7.0d0 ! nabla_y n first order
    density(2, 4, ipoint) = 8.0d0 ! nabla_z n first order
  end do

  res = derivative(id, &
                   num_grid_points, &
                   num_density_variables, &
                   vector_length, &
                   density)
  deallocate(density)

  ! compare with reference
  call assert((abs(-5.1509916226154067d0 - res) < 1.0d-6), &
              "derivatives do not match reference numbers")


  !-----------------------------------------------------------------------------
  ! now we will compute a particular partial derivative
  ! within order = 1
  ! we do this with a trick: we set the perturbed density for
  ! the density variable of interest to 1, and set other perturbed
  ! densities to 0

  allocate(density(vector_length, num_density_variables, num_grid_points))
  density = 0.0d0
  do ipoint = 1, num_grid_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 !         n zeroth order
    density(1, 2, ipoint) = 2.0d0 ! nabla_x n zeroth order
    density(1, 3, ipoint) = 3.0d0 ! nabla_y n zeroth order
    density(1, 4, ipoint) = 4.0d0 ! nabla_z n zeroth order
    density(2, 1, ipoint) = 0.0d0
    density(2, 2, ipoint) = 0.0d0
    density(2, 3, ipoint) = 1.0d0 ! we differentiate wrt this variable
    density(2, 4, ipoint) = 0.0d0
  end do

  res = derivative(id, &
                   num_grid_points, &
                   num_density_variables, &
                   vector_length, &
                   density)
  deallocate(density)

  ! compare with reference
  call assert((abs(-1.3470456737102541d-2 - res) < 1.0d-6), &
              "derivatives do not match reference numbers")


  !-----------------------------------------------------------------------------
  ! now we try 2nd order

  order = 2
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate(density(vector_length, num_density_variables, num_grid_points))
  density = 0.0d0
  do ipoint = 1, num_grid_points
    ! we use fantasy values here
    density(1, 1, ipoint) = 1.0d0 ! zeroth order
    density(1, 2, ipoint) = 2.0d0 ! zeroth order
    density(1, 3, ipoint) = 3.0d0 ! zeroth order
    density(1, 4, ipoint) = 4.0d0 ! zeroth order
    density(2, 1, ipoint) = 5.0d0 ! first order
    density(2, 2, ipoint) = 6.0d0 ! first order
    density(2, 3, ipoint) = 7.0d0 ! first order
    density(2, 4, ipoint) = 8.0d0 ! first order
    density(3, 1, ipoint) = 5.0d0 ! first order
    density(3, 2, ipoint) = 6.0d0 ! first order
    density(3, 3, ipoint) = 7.0d0 ! first order
    density(3, 4, ipoint) = 8.0d0 ! first order
    density(4, 1, ipoint) = 0.0d0 ! second order
    density(4, 2, ipoint) = 0.0d0 ! second order
    density(4, 3, ipoint) = 0.0d0 ! second order
    density(4, 4, ipoint) = 0.0d0 ! second order
  end do

  res = derivative(id, &
                   num_grid_points, &
                   num_density_variables, &
                   vector_length, &
                   density)
  deallocate(density)

  ! compare with reference
  call assert((abs(-9.4927931153398468d0 - res) < 1.0d-6), &
              "derivatives do not match reference numbers")


  !-----------------------------------------------------------------------------
  ! now we try 3nd order, contracted with perturbed densities

  order = 3
  ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
  call assert(ierr == 0, "xc_eval_setup failed")

  vector_length = 2**order
  allocate(density(vector_length, num_density_variables, num_grid_points))
  density = 0.0d0
  do ipoint = 1, num_grid_points
    ! we use fantasy values here
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

  res = derivative(id, &
                   num_grid_points, &
                   num_density_variables, &
                   vector_length, &
                   density)
  deallocate(density)

  ! compare with reference
  call assert((abs(47.091223089835331d0 - res) < 1.0d-6), &
              "derivatives do not match reference numbers")


  !-----------------------------------------------------------------------------
  ! we are done and can release the memory

  call xc_free_functional(id)

  print *, 'Kernel test has ended properly!'

contains

  real(8) function derivative(id, &
                              num_grid_points, &
                              num_density_variables, &
                              vector_length, &
                              density)
    ! computes the derivative and takes care of offsetting

    integer, intent(in) :: id
    integer, intent(in) :: num_grid_points
    integer, intent(in) :: num_density_variables
    integer, intent(in) :: vector_length
    real(8), intent(in) :: density(vector_length, num_density_variables, num_grid_points)

    real(8), allocatable :: input_array(:, :)
    real(8), allocatable :: output_array(:, :)

    integer :: ipoint

    allocate(input_array(num_density_variables*vector_length, num_grid_points))
    allocate(output_array(vector_length, num_grid_points))

    ! put the densities into the right places
    ! along the input array
    do ipoint = 1, num_grid_points
      input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                 density(:, 2, ipoint), &
                                 density(:, 3, ipoint), &
                                 density(:, 4, ipoint)/)
    end do

    call xc_eval(id, num_grid_points, input_array, output_array)

    ! The output_array holds a Taylor series expansion
    ! and we pick here one particular element out of this array.
    derivative = output_array(vector_length, 1)

    deallocate(input_array)
    deallocate(output_array)
  end function


  subroutine assert(predicate, error_message)

    logical, intent(in) :: predicate
    character(*), intent(in) :: error_message

    if (.not. predicate) then
      print *, 'ERROR:', error_message
      stop 1
    endif

  end subroutine

end program
