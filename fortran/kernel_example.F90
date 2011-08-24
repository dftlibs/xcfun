program xc_example

!  this example contains calls to all f90 interface routines
!  that are needed to "talk to" the xcfun library
!  and demonstrates how to use them

   use xcfun

   implicit none

   character(1000)      :: text
   integer              :: id, order, npoints, res, ideriv, i
   real(8), allocatable :: groundstate_density(:, :), output(:, :), density2(:,:,:), density3(:,:,:,:)
   real(8)              :: result_derv(4), result_derv_reference(4)

!  print some info and copyright about the library
!  please always include this info in your code
   call xcfun_splash(text)
   print *, text(1:len_trim(text))

!  create a new functional
!  we need this for interacting with the library
   id = xc_new_functional()

!  in this example we use PBE
   print *, 'Setting up PBE'
   call xc_set(id, XC_PBEX, 1.0d0)
   call xc_set(id, XC_PBEC, 1.0d0)
!  We will compute the kernel for an unpolarized system using total density and the
!  gradient components as the variables. These are linear in the density matrix, which
!  helps the code using the results from xcfun.

!  First we just compute the energy, i.e. the 0-th order integrand.
!  We have one gridpoint, and four variables, density N and gradient
!  components NX NY NZ.
   order = 0
   npoints = 1
   res = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
   if (res.ne.0) then
      print *,'xc_eval_setup failed with error ',res
      stop
   endif

!  Allocate space for the input
   allocate(groundstate_density(4,1))
!  and the output (one number, the XC energy density)
   allocate(output(1,1))
!  fill in some density
   groundstate_density(1,1) = 1.0 ! density
   groundstate_density(2,1) = 2.0 ! g_x
   groundstate_density(3,1) = 3.0 ! g_y
   groundstate_density(4,1) = 4.0 ! g_z
!  Ok, now compute the xc energy density at this point
   call xc_eval(id,npoints,groundstate_density,output)
   print *,'The XC energy density is',output(1,1)

   deallocate(groundstate_density)
   deallocate(output)

!  Now let's compute the first derivatives ('potential')
!  First set up xcfun for first derivatives
   order = 1
   res = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
   if (res.ne.0) then
      print *,'xc_eval_setup failed with error ',res
      stop
   endif

!  This time we need a two dimensional (+1 gridpoint dimension)
!  input, as follows
   allocate(density2(2,4,1))
!  and the output (two numbers in this "gridpoint", one energy, one derivative)
   allocate(output(2,1))
!  fill in the ground state density
   density2(1,1,1) = 1.0 ! density
   density2(1,2,1) = 2.0 ! g_x
   density2(1,3,1) = 3.0 ! g_y
   density2(1,4,1) = 4.0 ! g_z



   result_derv_reference(1) = -0.97182658347124051d0
   result_derv_reference(2) = -8.98025300594966838d-3
   result_derv_reference(3) = -1.34703795089245043d-2
   result_derv_reference(4) = -1.79605060118993368d-2

!  here is a tricky part, we need a loop and we need to put
!  in a '1' where we want the derivative
!  we loop over all variables we want a derivative with respect to
   do ideriv = 1, 4
      density2(2, :,      1) = 0.0d0
      density2(2, ideriv, 1) = 1.0d0
!     It is possible to put in other numbers than 1 and 0 in the density2 array
!     this is where perturbed densities go for automatic contraction

!     call xcfun to evaluate derivative with respect to the variable ideriv
      call xc_eval(id, npoints, reshape(density2, (/8, 1/)), output)
      print *, 'dE/dx_i for i =', ideriv, 'is', output(2, 1)

!     test against reference numbers
      result_derv(ideriv) = output(2, 1)
      if (abs(result_derv(ideriv) - result_derv_reference(ideriv)) > 1.0d-6) then
         print *, 'error: derivatives do not match reference numbers'
         stop 1
      end if
   end do

   deallocate(output)
!  Now second derivative of the "potential", contracted with one perturbed density.
!  Hopefully the strange input/output format will start to make sense
   order = 2
   res = xc_eval_setup(id, XC_N_NX_NY_NZ, XC_CONTRACTED, order)
   if (res.ne.0) then
      print *,'xc_eval_setup failed with error ',res
      stop
   endif
!  This time we need a three dimensional (+1 gridpoint dimension)
!  input, as follows
   allocate(density3(2,2,4,1))
!  and the output (four numbers in this "gridpoint", one energy, two first derivatives, one second derivative)
   allocate(output(4,1))
!  fill in the ground state density
   density3(1,1,1,1) = 1.0 ! density
   density3(1,1,2,1) = 2.0 ! g_x
   density3(1,1,3,1) = 3.0 ! g_y
   density3(1,1,4,1) = 4.0 ! g_z
!  fill in the perturbed density (a trial vector probably)
   density3(2,1,1,1) = 5.0 ! density
   density3(2,1,2,1) = 6.0 ! g_x
   density3(2,1,3,1) = 7.0 ! g_y
   density3(2,1,4,1) = 8.0 ! g_z

   result_derv_reference(1) = -2.0795504461938754d0
   result_derv_reference(2) =  2.14893341430147031d-2
   result_derv_reference(3) =  4.12142542204717230d-2
   result_derv_reference(4) =  6.09391742979287290d-2

!  now we have to put in ones and zeroes again to generate derivatives (which then generate matrix elements in ADF)
!  (1, 1, :) ground state values
!  (2, 1, :) perturbed density
!  (1, 2, :) ones and zeros
!  (2, 2, :) not used
!  in general there will be 2^N elements, 1 ground state, 2^N-2 perturbed and one unused.
   do ideriv = 1, 4
      density3(1, 2, :,      1) = 0.0d0
      density3(1, 2, ideriv, 1) = 1.0d0

!     call xcfun to evaluate derivative with respect to the variable ideriv
      call xc_eval(id, npoints, reshape(density3, (/16, 1/)), output)
      print *, 'd^2 E/dx_i dD_pert for i =', ideriv, 'is', output(4, 1)

!     test against reference numbers
      result_derv(ideriv) = output(4, 1)
      if (abs(result_derv(ideriv) - result_derv_reference(ideriv)) > 1.0d-6) then
         print *, 'error: derivatives do not match reference numbers'
         stop 1
      end if
   end do




   ! Note: the computed second derivatives are now already contracted with the
   ! perturbed density! This works to any order.
   ! In this example we put in once and zeros, which you need to do to compute matrix
   ! elements. You can also put in a perturbed density instead of ones and zeros,
   ! then you generate a response contribution d^2E/dD1dD2. Typically this is what
   ! you want to do to high order, because by the 2N+1 rule you don't have to compute
   ! matrix elements to very high order. Then you need only one call of xc_eval,
   ! not one for each variable. This is very efficient, much better than computing
   ! all partial derivatives. For second derivative it might be better to compute
   ! all partial derivatives and reuse them in the response solved.

   ! Note2: In this example I used the variables that are linear in the density
   ! matrix. This makes life easier because you can construct them trivially from
   ! perturbed density matrices. If you use i.e. the square norm of the density
   ! gradient things get more complicated.
   ! Note3: You can extend this example trivially to alpha/beta densities
   ! By specifying XC_A_B_AX_AY_AZ_BX_BY_BZ instead of  XC_N_NX_NY_NZ
   ! Then the number 4 above will be replaced by 8.

end program
