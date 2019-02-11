#include "../src/cpp_macros.h"
!> \example dielectric_Xd.f90
!>
!> Example showing how to include a dielectric object. Warning: the
!> functionality is not fully ready
program dielectric_test
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_eps

  ! The dielectric constant used in this example
  double precision, parameter :: epsilon_high = 10.0_dp

  type(af_t)        :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter
  real(dp)           :: residu(2)
  character(len=100) :: fname
  type(mg_t)        :: mg
  integer            :: count_rate, t_start, t_end

  print *, "****************************************"
  print *, "Warning: functionality demonstrated here is not fully ready"
  print *, "For large epsilon, convergence will probably be slow"
  print *, "****************************************"
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "eps", ix=i_eps)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)])

  call af_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call af_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call af_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  ! Average epsilon on coarse grids. In the future, it could be better to define
  ! epsilon on cell faces, and to perform this restriction in a matrix fashion:
  ! A_coarse = M_restrict * A_fine * M_prolong (A = matrix operator, M = matrix)
  call af_restrict_tree(tree, i_eps)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call af_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%i_eps        = i_eps       ! Variable for epsilon coefficient
  mg%sides_bc     => sides_bc

  ! Automatically detect the right methods
  mg%box_op       => mg_auto_op
  mg%box_gsrb     => mg_auto_gsrb
  mg%box_corr     => mg_auto_corr
  mg%box_stencil  => mg_box_lpld_stencil

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg_fas_fmg
  call mg_init(tree, mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)

     ! Determine the minimum and maximum residual and error
     call af_tree_min_cc(tree, i_tmp, residu(1))
     call af_tree_max_cc(tree, i_tmp, residu(2))
     write(*,"(I8,Es14.5)") mg_iter, maxval(abs(residu))

     write(fname, "(A,I0)") "dielectric_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call af_destroy(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    real(dp)                 :: eps_min, eps_max

    eps_min = minval(box%cc(DTIMES(:), i_eps))
    eps_max = maxval(box%cc(DTIMES(:), i_eps))

    if ((box%lvl < 5 .and. eps_max > eps_min) .or. box%lvl < 2) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                      :: IJK, nc
    real(dp)                     :: rr(NDIM)
    real(dp)                     :: ellips_fac(NDIM)

    nc = box%n_cell

    ! Create ellipsoidal shape
    ellips_fac(2:) = 3.0_dp
    ellips_fac(1)  = 1.0_dp

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])

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
    use m_af_ghostcell
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (af_neighb_lowx)
       call af_bc_dirichlet_zero(box, nb, iv, bc_type)
    case (af_neighb_highx)
       bc_type = af_bc_dirichlet
#if NDIM == 2
       box%cc(nc+1, 1:nc, iv) = 1.0_dp
#elif NDIM == 3
       box%cc(nc+1, 1:nc, 1:nc, iv) = 1.0_dp
#endif
    case default
       call af_bc_neumann_zero(box, nb, iv, bc_type)
    end select
  end subroutine sides_bc

end program dielectric_test
