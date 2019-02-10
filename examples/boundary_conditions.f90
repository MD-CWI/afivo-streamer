#include "../src/cpp_macros.h"
!> \example boundary_conditions_Xd.f90
!>
!> Example showing how to use different types of boundary conditions
program boundary_conditions_Xd
  use m_af_all

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 20
  integer            :: i_phi

  type(af_t)         :: tree
  integer            :: iter
  character(len=100) :: fname

  print *, "Running boundary_conditions_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)])

  call af_print_info(tree)

  do iter = 1, n_iterations
     if (iter == 1) then
        ! Set initial conditions on first iteration
        call af_loop_box(tree, set_phi_zero)
     else
        ! Replace phi by average over neighbors
        call af_loop_box(tree, average_phi)
     end if

     call af_gc_tree(tree, i_phi, af_gc_interp, boundary_method)

     write(fname, "(A,I0)") "boundary_conditions_" // DIMNAME // "_", iter
     call af_write_vtk(tree, trim(fname), dir="output")
  end do

contains

  ! This routine sets the initial conditions for each box
  subroutine set_phi_zero(box)
    type(box_t), intent(inout) :: box

    box%cc(DTIMES(:), i_phi) = 0.0_dp
  end subroutine set_phi_zero

  subroutine average_phi(box)
    type(box_t), intent(inout) :: box
    integer                      :: IJK, nc
    real(dp)                     :: tmp(DTIMES(box%n_cell))

    nc = box%n_cell

    do KJI_DO(1,nc)
#if NDIM == 2
       tmp(i, j) = 0.25_dp * ( &
            box%cc(i+1, j, i_phi) + &
            box%cc(i-1, j, i_phi) + &
            box%cc(i, j+1, i_phi) + &
            box%cc(i, j-1, i_phi))
#elif NDIM == 3
       tmp(i, j, k) = 1/6.0_dp * ( &
            box%cc(i+1, j, k, i_phi) + &
            box%cc(i-1, j, k, i_phi) + &
            box%cc(i, j+1, k, i_phi) + &
            box%cc(i, j-1, k, i_phi) + &
            box%cc(i, j, k+1, i_phi) + &
            box%cc(i, j, k-1, i_phi))
#endif
    end do; CLOSE_DO

    box%cc(DTIMES(1:nc), i_phi) = tmp(DTIMES(:))
  end subroutine average_phi

  !> [boundary_method]
  subroutine boundary_method(box, nb, iv, bc_type)
    type(box_t), intent(inout) :: box     ! Box to operate on
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc

    nc = box%n_cell

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)      ! Lower-x direction
       bc_type = af_bc_dirichlet
       box%cc(0, 1:nc, iv) = 1.0_dp
    case (af_neighb_highx)     ! Higher-x direction
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0.0_dp
    case (af_neighb_lowy)      ! Lower-y direction
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 0, iv) = 1.0_dp
    case (af_neighb_highy)     ! Higher-y direction
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, iv) = 0.0_dp
#elif NDIM == 3
    case (af_neighb_lowx)      ! Lower-x direction
       bc_type = af_bc_dirichlet
       box%cc(0, 1:nc, 1:nc, iv) = 1.0_dp
    case (af_neighb_highx)     ! Higher-x direction
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, 1:nc, iv) = 0.0_dp
    case (af_neighb_lowy)      ! Lower-y direction
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 0, 1:nc, iv) = 1.0_dp
    case (af_neighb_highy)     ! Higher-y direction
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, 1:nc, iv) = 0.0_dp
    case (af_neighb_lowz)      ! Lower-z direction
       bc_type = af_bc_neumann
       box%cc(1:nc, 1:nc, 0, iv) = 0.0_dp
    case (af_neighb_highz)     ! Higher-z direction
       bc_type = af_bc_neumann
       box%cc(1:nc, 1:nc, nc+1, iv) = 0.0_dp
#endif
    end select

  end subroutine boundary_method
  !> [boundary_method]

end program boundary_conditions_Xd
