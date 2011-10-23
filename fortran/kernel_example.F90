program xc_example

!  this example contains calls to all f90 interface routines
!  that are needed to "talk to" the xcfun library
!  and demonstrates how to use them
!  We will compute the kernel for an unpolarized system using total density and the
!  gradient components as the variables. These are linear in the density matrix, which
!  helps the code using the results from xcfun.

   use xcfun

   implicit none

   integer, parameter   :: NR_POINTS = 1
   integer, parameter   :: NR_VARIABLES = 4

   character(1000)      :: text
   integer              :: id, order, ierr, ideriv, ipoint
   integer              :: vector_length
   real(8)              :: res, res_reference(NR_VARIABLES)

   real(8), allocatable :: density(:, :, :)
   real(8), allocatable :: input_array(:, :)
   real(8), allocatable :: output_array(:, :)

!  print some info and copyright about the library
!  please always include this info in your code
   call xcfun_splash(text)
   print *, text(1:len_trim(text))

!  create a new functional
!  we need this for interacting with the library
   id = xc_new_functional()

!  in this example we use PBE
   print *, 'Setting up PBE'
   ierr = xc_set(id, 'pbe', 1.0d0)
   if (ierr.ne.0) then
      print *, 'Functional name not recognized, quitting.'
      stop
   endif

!  First we just compute the energy, i.e. the 0-th order integrand.
!  We have one gridpoint, and four variables, density N and gradient
!  components NX NY NZ.
   order = 0
   ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
   if (ierr /= 0) then
      print *, 'xc_eval_setup failed with error ', ierr
      stop 1
   endif

   vector_length = 2**order
!  allocate density and fill in some phantasy values
   allocate(density(vector_length, NR_VARIABLES, NR_POINTS))
   density = 0.0d0
   do ipoint = 1, NR_POINTS
      density(1, 1, ipoint) = 1.0d0 !        n
      density(1, 2, ipoint) = 2.0d0 !nabla_x n
      density(1, 3, ipoint) = 3.0d0 !nabla_y n
      density(1, 4, ipoint) = 4.0d0 !nabla_z n
   end do

   allocate(input_array(NR_VARIABLES*vector_length, NR_POINTS))
   allocate(output_array(vector_length, NR_POINTS))

   do ipoint = 1, NR_POINTS
      input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                 density(:, 2, ipoint), &
                                 density(:, 3, ipoint), &
                                 density(:, 4, ipoint)/)
   end do
   call xc_eval(id, NR_POINTS, input_array, output_array)
   res = output_array(vector_length, 1)
   print *, 'The XC energy density is', res

   deallocate(density)
   deallocate(input_array)
   deallocate(output_array)

!  reference values for self test
   res_reference(1) = -0.97182658347124051d0
   res_reference(2) = -8.98025300594966838d-3
   res_reference(3) = -1.34703795089245043d-2
   res_reference(4) = -1.79605060118993368d-2

!  now let's compute the first derivatives ('potential')
!  first set up xcfun for first derivatives
   order = 1
   ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
   if (ierr /= 0) then
      print *, 'xc_eval_setup failed with error ', ierr
      stop 1
   endif

   vector_length = 2**order
!  allocate density and fill in some phantasy values
   allocate(density(vector_length, NR_VARIABLES, NR_POINTS))
   density = 0.0d0
   do ipoint = 1, NR_POINTS
      density(1, 1, ipoint) = 1.0d0 !        n
      density(1, 2, ipoint) = 2.0d0 !nabla_x n
      density(1, 3, ipoint) = 3.0d0 !nabla_y n
      density(1, 4, ipoint) = 4.0d0 !nabla_z n
   end do

   allocate(input_array(NR_VARIABLES*vector_length, NR_POINTS))
   allocate(output_array(vector_length, NR_POINTS))

!  here is a tricky part, we need a loop and we need to put
!  in a '1' where we want the derivative
!  we loop over all variables we want a derivative with respect to
   do ideriv = 1, NR_VARIABLES
      do ipoint = 1, NR_POINTS
         density(2, :,      ipoint) = 0.0d0
         density(2, ideriv, ipoint) = 1.0d0
      end do
!     It is possible to put in other numbers than 1 and 0
!     this is where perturbed densities go for automatic contraction

      do ipoint = 1, NR_POINTS
         input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                    density(:, 2, ipoint), &
                                    density(:, 3, ipoint), &
                                    density(:, 4, ipoint)/)
      end do
      call xc_eval(id, NR_POINTS, input_array, output_array)
      res = output_array(vector_length, 1)
      print *, 'dE/dx_i for i =', ideriv, 'is', res

!     test against reference numbers
      if (abs(res - res_reference(ideriv)) > 1.0d-6) then
         print *, 'error: derivatives do not match reference numbers'
         stop 1
      end if
   end do

   deallocate(density)
   deallocate(input_array)
   deallocate(output_array)

!  reference values for self test
   res_reference(1) = -2.0795504461938754d0
   res_reference(2) =  2.14893341430147031d-2
   res_reference(3) =  4.12142542204717230d-2
   res_reference(4) =  6.09391742979287290d-2

!  now second derivative of the "potential", contracted with one perturbed density.
!  hopefully the strange input/output_array format will start to make sense
   order = 2
   ierr = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
   if (ierr /= 0) then
      print *, 'xc_eval_setup failed with error ', ierr
      stop 1
   endif

   vector_length = 2**order
!  allocate density and fill in some phantasy values
   allocate(density(vector_length, NR_VARIABLES, NR_POINTS))
   density = 0.0d0
   do ipoint = 1, NR_POINTS
      density(1, 1, ipoint) = 1.0d0
      density(1, 2, ipoint) = 2.0d0
      density(1, 3, ipoint) = 3.0d0
      density(1, 4, ipoint) = 4.0d0
      density(2, 1, ipoint) = 5.0d0
      density(2, 2, ipoint) = 6.0d0
      density(2, 3, ipoint) = 7.0d0
      density(2, 4, ipoint) = 8.0d0
   end do

   allocate(input_array(NR_VARIABLES*vector_length, NR_POINTS))
   allocate(output_array(vector_length, NR_POINTS))

!  now we have to put in ones and zeroes again to generate derivatives (which then generate matrix elements in ADF)
!  density(1, :, :)  unperturbed values
!  density(2, :, :)  perturbed values
!  density(3, :, :)  ones and zeros
!  density(4, :, :)  not used
!  in general there will be 2^order elements, 1 ground state, 2^order-2 perturbed and one unused.
   do ideriv = 1, NR_VARIABLES
      do ipoint = 1, NR_POINTS
         density(3, :,      ipoint) = 0.0d0
         density(3, ideriv, ipoint) = 1.0d0
      end do
!     It is possible to put in other numbers than 1 and 0 in the density2 array
!     this is where perturbed densities go for automatic contraction

      do ipoint = 1, NR_POINTS
         input_array(:, ipoint) = (/density(:, 1, ipoint), &
                                    density(:, 2, ipoint), &
                                    density(:, 3, ipoint), &
                                    density(:, 4, ipoint)/)
      end do
      call xc_eval(id, NR_POINTS, input_array, output_array)
      res = output_array(vector_length, 1)
      print *, 'd^2 E/dx_i dD_pert for i =', ideriv, 'is', res

!     test against reference numbers
      if (abs(res - res_reference(ideriv)) > 1.0d-6) then
         print *, 'error: derivatives do not match reference numbers'
         stop 1
      end if
   end do

   deallocate(density)
   deallocate(input_array)
   deallocate(output_array)

!  Note: the computed second derivatives are now already contracted with the
!  perturbed density! This works to any order.
!  In this example we put in ones and zeros, which you need to do to compute matrix
!  elements. You can also put in a perturbed density instead of ones and zeros,
!  then you generate a response contribution d^2E/dD1dD2. Typically this is what
!  you want to do to high order, because by the 2N+1 rule you don't have to compute
!  matrix elements to very high order. Then you need only one call of xc_eval,
!  not one for each variable. This is very efficient, much better than computing
!  all partial derivatives. For second derivatives it might be better to compute
!  all partial derivatives and reuse them in the response solver.

!  Note2: In this example I used the variables that are linear in the density
!  matrix. This makes life easier because you can construct them trivially from
!  perturbed density matrices. If you use i.e. the square norm of the density
!  gradient things get more complicated.

!  Note3: You can extend this example trivially to alpha/beta densities
!  By specifying XC_A_B_AX_AY_AZ_BX_BY_BZ instead of  XC_N_NX_NY_NZ
!  Then the NR_VARIABLES becomes 8.

end program
