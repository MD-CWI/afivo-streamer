#include "../src/cpp_macros_$Dd.h"
!> \example boundary_conditions_$Dd.f90
!>
!> Example showing how to use different types of boundary conditions
program boundary_conditions_$Dd
  use m_a$D_all

  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 20
  integer, parameter :: i_phi = 1

  type(a$D_t)        :: tree
  integer            :: iter
  integer            :: ix_list($D, n_boxes_base)
  real(dp)           :: dr
  character(len=100) :: fname

  print *, "Running boundary_conditions_$Dd"
  print *, "Number of threads", af_get_max_threads()

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  ! Initialize tree
  call a$D_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       1, &            ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       cc_names=["phi"]) ! Variable names

  ! Set up geometry. This creates a single box, for which we need boundary
  ! conditions on all sides
  ix_list(:, 1) = [DTIMES(1)]         ! Set index of box 1

  ! Create the base mesh, using the box indices and their neighbor information
  call a$D_set_base(tree, 1, ix_list)
  call a$D_print_info(tree)

  do iter = 1, n_iterations
     if (iter == 1) then
        ! Set initial conditions on first iteration
        call a$D_loop_box(tree, set_phi_zero)
     else
        ! Replace phi by average over neighbors
        call a$D_loop_box(tree, average_phi)
     end if

     call a$D_gc_tree(tree, i_phi, a$D_gc_interp, boundary_method)

     write(fname, "(A,I0)") "boundary_conditions_$Dd_", iter
     call a$D_write_vtk(tree, trim(fname), dir="output")
  end do

contains

  ! This routine sets the initial conditions for each box
  subroutine set_phi_zero(box)
    type(box$D_t), intent(inout) :: box

    box%cc(DTIMES(:), i_phi) = 0.0_dp
  end subroutine set_phi_zero

  subroutine average_phi(box)
    type(box$D_t), intent(inout) :: box
    integer                      :: IJK, nc
    real(dp)                     :: tmp(DTIMES(box%n_cell))

    nc = box%n_cell

    do KJI_DO(1,nc)
#if $D == 2
       tmp(i, j) = 0.25_dp * ( &
            box%cc(i+1, j, i_phi) + &
            box%cc(i-1, j, i_phi) + &
            box%cc(i, j+1, i_phi) + &
            box%cc(i, j-1, i_phi))
#elif $D == 3
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
    type(box$D_t), intent(inout) :: box     ! Box to operate on
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc

    nc = box%n_cell

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
#if $D == 2
    case (a$D_neighb_lowx)      ! Lower-x direction
       bc_type = af_bc_dirichlet
       box%cc(0, 1:nc, iv) = 1.0_dp
    case (a$D_neighb_highx)     ! Higher-x direction
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0.0_dp
    case (a$D_neighb_lowy)      ! Lower-y direction
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 0, iv) = 1.0_dp
    case (a$D_neighb_highy)     ! Higher-y direction
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, iv) = 0.0_dp
#elif $D == 3
    case (a$D_neighb_lowx)      ! Lower-x direction
       bc_type = af_bc_dirichlet
       box%cc(0, 1:nc, 1:nc, iv) = 1.0_dp
    case (a$D_neighb_highx)     ! Higher-x direction
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, 1:nc, iv) = 0.0_dp
    case (a$D_neighb_lowy)      ! Lower-y direction
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 0, 1:nc, iv) = 1.0_dp
    case (a$D_neighb_highy)     ! Higher-y direction
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, 1:nc, iv) = 0.0_dp
    case (a$D_neighb_lowz)      ! Lower-z direction
       bc_type = af_bc_neumann
       box%cc(1:nc, 1:nc, 0, iv) = 0.0_dp
    case (a$D_neighb_highz)     ! Higher-z direction
       bc_type = af_bc_neumann
       box%cc(1:nc, 1:nc, nc+1, iv) = 0.0_dp
#endif
    end select

  end subroutine boundary_method
  !> [boundary_method]

end program boundary_conditions_$Dd
