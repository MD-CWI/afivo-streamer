#include "../src/cpp_macros.h"
!> \example boundary_conditions.f90
!>
!> Example showing how to use different types of boundary conditions
program boundary_conditions
  use m_af_all

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 20
  integer            :: i_phi, i_rho

  type(af_t)         :: tree
  integer            :: iter
  character(len=100) :: fname

  print *, "Running boundary_conditions_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rho", ix=i_rho)
  call af_set_cc_methods(tree, i_phi, boundary_method)
  call af_set_cc_methods(tree, i_rho, bc_custom=custom_boundary_method)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)])

  call af_print_info(tree)

  do iter = 1, n_iterations
     if (iter == 1) then
        ! Set initial conditions on first iteration
        call af_loop_box(tree, set_initial_cond)
     else
        ! Replace phi by average over neighbors
        call af_loop_box(tree, average_vars)
     end if

     call af_gc_tree(tree, [i_phi, i_rho])

     write(fname, "(A,I0)") "output/boundary_conditions_" // DIMNAME // "_", iter
     call af_write_silo(tree, trim(fname))
  end do

contains

  ! This routine sets the initial conditions for each box
  subroutine set_initial_cond(box)
    type(box_t), intent(inout) :: box

    box%cc(DTIMES(:), i_phi) = 0.0_dp
    box%cc(DTIMES(:), i_rho) = 1.0_dp
  end subroutine set_initial_cond

  subroutine average_vars(box)
    type(box_t), intent(inout) :: box
    integer                      :: IJK, nc, ivs(2)
    real(dp)                     :: tmp(DTIMES(box%n_cell), 2)

    nc = box%n_cell
    ivs = [i_phi, i_rho]

    do KJI_DO(1,nc)
#if NDIM == 1
       tmp(i, :) = 0.5_dp * ( &
            box%cc(i+1, ivs) + &
            box%cc(i-1, ivs))
#elif NDIM == 2
       tmp(i, j, :) = 0.25_dp * ( &
            box%cc(i+1, j, ivs) + &
            box%cc(i-1, j, ivs) + &
            box%cc(i, j+1, ivs) + &
            box%cc(i, j-1, ivs))
#elif NDIM == 3
       tmp(i, j, k, :) = 1/6.0_dp * ( &
            box%cc(i+1, j, k, ivs) + &
            box%cc(i-1, j, k, ivs) + &
            box%cc(i, j+1, k, ivs) + &
            box%cc(i, j-1, k, ivs) + &
            box%cc(i, j, k+1, ivs) + &
            box%cc(i, j, k-1, ivs))
#endif
    end do; CLOSE_DO

    ! Average new and old value
    box%cc(DTIMES(1:nc), ivs) = 0.5_dp * (&
         tmp(DTIMES(:), :) + box%cc(DTIMES(1:nc), ivs))
  end subroutine average_vars

  !> [boundary_method]
  subroutine boundary_method(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: nc

    nc = box%n_cell

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
    case (af_neighb_lowx)      ! Lower-x direction
       bc_type = af_bc_dirichlet
       bc_val = 1.0_dp
    case default
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end select

  end subroutine boundary_method
  !> [boundary_method]

  !> With this method we can set ghost cells manually
  subroutine custom_boundary_method(box, nb, iv, n_gc, cc)
    type(box_t), intent(inout) :: box     !< Box that needs b.c.
    integer, intent(in)     :: nb      !< Direction
    integer, intent(in)     :: iv      !< Index of variable
    integer, intent(in)     :: n_gc !< Number of ghost cells to fill
    !> If n_gc > 1, fill ghost cell values in this array instead of box%cc
    real(dp), intent(inout), optional :: &
         cc(DTIMES(1-n_gc:box%n_cell+n_gc))

    integer :: lo(NDIM), hi(NDIM)

    if (n_gc /= 1) error stop "not implemented"

    ! Get ghost cell index range
    call af_get_index_bc_outside(nb, box%n_cell, 1, lo, hi)

    ! Set all ghost cells to a scalar value
    box%cc(DSLICE(lo,hi), iv) = 0.0_dp
  end subroutine custom_boundary_method

end program
