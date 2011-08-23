program xc_example

!  this example contains all f90 interface routines
!  that are needed to "talk to" the xcfun library
!  and demonstrates how to use them

   use xcfun

   implicit none

   character(1000)      :: text
   integer              :: id, order, ilen, olen, res
   integer              :: i, k, ipoint, nr_points, block_length, max_block_length
   real(8)              :: derivative_nn_ab, derivative_ss_ab
   real(8)              :: derivative_nn_rs, derivative_ss_rs
   real(8), allocatable :: density_variables(:, :), derivatives(:, :)

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

!  print currently set parameters
!  discover all settings
!   print *, 'Now defined (nonzero weights and parameters, including defaults):'
!   do i = 1, XC_NR_PARAMS
!      call xc_param_name(i, text)
!      if (dabs(xc_get_param(id, i)) > tiny(0.0d0)) then
!         print *, text(4:len_trim(text)), xc_get_param(id, i)
!      end if
!   end do

!  functional type
!   print *, 'Functional type:', xc_get_type(id), ' (0 LDA, 1 GGA, 2 MGGA)'

!  let's get derivatives to second order
   order = 2
   res = xc_eval_setup(id, XC_A_B_GAA_GAB_GBB, XC_PARTIAL_DERIVATIVES, order)
   if (res.eq.0) then
      print *,'xc_eval_setup ok ',res
   else
      print *,'xc_eval_setup failed with error ',res
   endif

#if 0
   print *, 'Order:', order
   if (.not.xc_try_order(id,order)) then
      print *, 'Could not set order ',order
      stop
   endif

   if (.not.xc_try_vars(id,XC_A_B_GAA_GAB_GBB)) then
      print *, 'Could not set order ',order
      stop
   endif

   olen = xc_output_length(id, order)
   print *, 'Length of output (how many derivatives):', olen

!  set alpha/beta density variable mode
   call xc_set_mode(id, XC_VARS_AB)

!  here we decide to evaluate derivatives at 25 points
!  by batches of 10 (or less; for the last batch)
   nr_points        = 25
   max_block_length = 10

!  allocate space for density variables
   allocate(density_variables(ilen, max_block_length))
!  allocate space for functional derivatives
   allocate(derivatives(olen, max_block_length))

   do ipoint = 1, nr_points, max_block_length

      block_length = min(max_block_length, nr_points - ipoint + 1)

!     evaluate density variables in point(s)
!     (here we cheat and get some random numbers)
      do k = 1, block_length
         do i = 1, ilen
            call random_number(density_variables(i, k))
         end do
      end do

!     compute the functional and its derivatives up to order
!     we use ilen and olen as the leading dimension (pitch) values
!     although with a single point these are not important
      call xc_eval(id, order, block_length, density_variables, derivatives)

!     in this example we are interested in say
!     (d/dna)^1 (d/dnb)^0 (d/dgaa)^1 (d/dgab)^0 (d/dgab)^0 e_xc
      print *, 'An example second derivative'
      print *, '(d/dna)^1 (d/dnb)^0 (d/dgaa)^1 (d/dgab)^0 (d/dgab)^0 e_xc'
      do k = 1, block_length
         print *, 'k, exc_10100:', k, derivatives(XC_D10100, k)
      end do
      print *, 'Alternatively by xc_index (a bit slower for large grids):'
      do k = 1, block_length
         print *, 'k, exc_10100:', k, derivatives(xc_index(id, (/1, 0, 1, 0, 0/)), k)
      end do

   end do

!  deallocate arrays
   deallocate(density_variables)
   deallocate(derivatives)

!  finally let's compare the variable modes XC_VARS_AB (alpha and beta densities)
!  and the complementary XC_VARS_NS (total density and spin density mode)
!  at closed-shell reference for LDA

   id = xc_new_functional()
   call xc_set_param(id, XC_SLATERX, 1.0d0)
   call xc_set_param(id, XC_VWN5C,   1.0d0)
   ilen = xc_input_length(id)
   olen = xc_output_length(id, order)

   allocate(density_variables(ilen, 1))
   allocate(derivatives(olen, 1))

!  we switch to total/spin density mode
   call xc_set_mode(id, XC_VARS_NS)

   density_variables(1, 1) = 0.4d0
   density_variables(2, 1) = 0.0d0

   call xc_eval(id, order, 1, density_variables, derivatives)
   derivative_nn_rs = derivatives(XC_D20, 1)
   derivative_ss_rs = derivatives(XC_D02, 1)

!  we switch back to alpha/beta density mode
   call xc_set_mode(id, XC_VARS_AB)

   density_variables(1, 1) = 0.2d0
   density_variables(2, 1) = 0.2d0

   call xc_eval(id, order, 1, density_variables, derivatives)
   derivative_nn_ab = 0.5d0*derivatives(XC_D20, 1) + 0.5d0*derivatives(XC_D11, 1)
   derivative_ss_ab = 0.5d0*derivatives(XC_D20, 1) - 0.5d0*derivatives(XC_D11, 1)

   print *, 'Comparison of derivatives via XC_VARS_NS (left) and XC_VARS_AB (right):'
   print *, '(d/dn)^2 (d/ds)^0 e_xc:', derivative_nn_rs, derivative_nn_ab
   print *, '(d/dn)^0 (d/ds)^2 e_xc:', derivative_ss_rs, derivative_ss_ab

   deallocate(density_variables)
   deallocate(derivatives)
#endif
end program
