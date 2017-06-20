#include "../src/cpp_macros_$Dd.h"
!> \example dielectric_$Dd.f90
!>
!> Example showing how to include a dielectric object. Warning: the
!> functionality is not fully ready
program dielectric_test
  use m_a$D_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 10
  integer, parameter :: n_var_cell = 4
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_tmp = 3
  integer, parameter :: i_eps = 4

  ! The dielectric constant used in this example
  double precision, parameter :: epsilon_high = 10.0_dp

  type(a$D_t)        :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter
  integer            :: ix_list($D, n_boxes_base)
  real(dp)           :: dr, residu(2)
  character(len=100) :: fname
  type(mg$D_t)        :: mg
  integer            :: count_rate, t_start, t_end

  print *, "****************************************"
  print *, "Warning: functionality demonstrated here is not fully ready"
  print *, "For large epsilon, convergence will probably be slow"
  print *, "****************************************"
  print *, "Number of threads", af_get_max_threads()

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  ! Initialize tree
  call a$D_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "tmp", "eps"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = 1             ! Set index of box 1

  ! Create the base mesh, using the box indices and their neighbor information
  call a$D_set_base(tree, 1, ix_list)
  call a$D_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call a$D_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a$D_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  ! Average epsilon on coarse grids. In the future, it could be better to define
  ! epsilon on cell faces, and to perform this restriction in a matrix fashion:
  ! A_coarse = M_restrict * A_fine * M_prolong (A = matrix operator, M = matrix)
  call a$D_restrict_tree(tree, i_eps)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call a$D_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%i_eps        = i_eps       ! Variable for epsilon coefficient
  mg%sides_bc     => sides_bc

  ! Automatically detect the right methods
  mg%box_op       => mg$D_auto_op
  mg%box_gsrb     => mg$D_auto_gsrb
  mg%box_corr     => mg$D_auto_corr

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg$D_fas_fmg
  call mg$D_init_mg(mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg$D_fas_fmg(tree, mg, .true., mg_iter>1)

     ! Determine the minimum and maximum residual and error
     call a$D_tree_min_cc(tree, i_tmp, residu(1))
     call a$D_tree_max_cc(tree, i_tmp, residu(2))
     write(*,"(I8,Es14.5)") mg_iter, maxval(abs(residu))

     write(fname, "(A,I0)") "dielectric_$Dd_", mg_iter
     call a$D_write_silo(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a$D_destroy(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box$D_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    real(dp)                 :: eps_min, eps_max

    eps_min = minval(box%cc(DTIMES(:), i_eps))
    eps_max = maxval(box%cc(DTIMES(:), i_eps))

    if ((box%lvl < 7 .and. eps_max > eps_min) .or. box%lvl < 3) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box$D_t), intent(inout) :: box
    integer                      :: IJK, nc
    real(dp)                     :: rr($D)
    real(dp)                     :: ellips_fac($D)

    nc = box%n_cell
    dr = box%dr

    ! Create ellipsoidal shape
    ellips_fac(2:) = 3.0_dp
    ellips_fac(1)  = 1.0_dp

    do KJI_DO(0,nc+1)
       rr = a$D_r_cc(box, [IJK])

       ! Change epsilon in part of the domain
       if (norm2((rr - 0.5_dp) * ellips_fac) < 0.25_dp) then
          box%cc(IJK, i_eps) = epsilon_high
       else
          box%cc(IJK, i_eps) = 1.0_dp
       end if

       box%cc(IJK, i_rhs) = 0.0d0
       box%cc(IJK, i_phi) = 0.0d0
    end do; CLOSE_DO

  end subroutine set_init_cond

  subroutine sides_bc(box, nb, iv, bc_type)
    use m_a$D_ghostcell
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (a$D_neighb_lowx)
       call a$D_bc_dirichlet_zero(box, nb, iv, bc_type)
    case (a$D_neighb_highx)
       bc_type = af_bc_dirichlet
#if $D == 2
       box%cc(nc+1, 1:nc, iv) = 1.0_dp
#elif $D == 3
       box%cc(nc+1, 1:nc, 1:nc, iv) = 1.0_dp
#endif
    case default
       call a$D_bc_neumann_zero(box, nb, iv, bc_type)
    end select
  end subroutine sides_bc

end program dielectric_test
