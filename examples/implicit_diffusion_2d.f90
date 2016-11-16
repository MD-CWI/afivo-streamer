!> \example drift_diffusion_2d.f90
!> A drift-diffusion example
!> @TODO: document this
program drift_diffusion_2d
  use m_a2_types
  use m_a2_core
  use m_a2_ghostcell
  use m_a2_utils
  use m_a2_output
  use m_a2_restrict
  use m_a2_multigrid

  implicit none

  integer, parameter :: box_size    = 8
  integer, parameter :: i_phi       = 1
  integer, parameter :: i_rhs       = 2
  integer, parameter :: i_tmp       = 3
  integer, parameter :: i_err       = 4

  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  real(dp), parameter :: dr = domain_len / box_size
  real(dp), parameter :: diffusion_coeff = 0.1_dp

  type(a2_t)         :: tree
  type(mg2_t)        :: mg
  type(ref_info_t)   :: ref_info
  integer            :: id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer            :: time_steps, output_cnt
  real(dp)           :: dt, time, end_time
  character(len=100) :: fname

  print *, "Running drift_diffusion_2d"
  print *, "Number of threads", af_get_max_threads()

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=4, n_var_face=1, dr=dr, &
       cc_names=["phi", "rhs", "tmp", "err"])

  ! Set up geometry
  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = 1     ! Periodic domain

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)
  call a2_print_info(tree)

  ! Set up the initial conditions
  do
     call a2_loop_box(tree, set_initial_condition)
     call a2_gc_tree(tree, i_phi, a2_gc_interp, a2_bc_dirichlet_zero)
     call a2_adjust_refinement(tree, refine_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  mg%i_phi    = i_phi                 ! Solution variable
  mg%i_rhs    = i_rhs                 ! Right-hand side variable
  mg%i_tmp    = i_tmp                 ! Variable for temporary space
  mg%sides_bc => a2_bc_dirichlet_zero ! Method for boundary conditions

  ! The methods defined below implement a backward Euler method for the heat
  ! equation, by changing the elliptic operator for the multigrid procedure.
  mg%box_op   => box_op_diff
  mg%box_gsrb => box_gsrb_diff

  ! This routine does not initialize the multigrid fields boxes%i_phi,
  ! boxes%i_rhs and boxes%i_tmp. These fileds will be initialized at the
  ! first call of mg2_fas_fmg
  call mg2_init_mg(mg)

  output_cnt = 0
  time       = 0
  end_time   = 2.0_dp
  time_steps = 0
  time       = 0
  dt         = 0.1_dp

  ! Starting simulation
  do
     time_steps = time_steps + 1

     output_cnt = output_cnt + 1
     write(fname, "(A,I0)") "implicit_diffusion_2d_", output_cnt
     call a2_loop_box_arg(tree, set_error, [time])
     call a2_write_vtk(tree, trim(fname), output_cnt, time, dir="output")

     if (time > end_time) exit

     call a2_loop_box(tree, set_rhs)
     call mg2_fas_fmg(tree, mg, set_residual=.true., have_guess=.true.)
     time = time + dt
  end do

contains

  ! Return the refinement flag for boxes(id)
  subroutine refine_routine(boxes, id, refine_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: refine_flag

    if (boxes(id)%dr > 5e-3_dp * domain_len) refine_flag = af_do_ref
  end subroutine refine_routine

  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_phi) = solution(xy, 0.0_dp)
       end do
    end do
  end subroutine set_initial_condition

  subroutine set_error(box, time)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: time(:)
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = &
               box%cc(i, j, i_phi) - solution(xy, time(1))
       end do
    end do
  end subroutine set_error

  function solution(xy, t) result(sol)
    real(dp), intent(in) :: xy(2), t
    real(dp)             :: sol
    integer, parameter   :: cxy(2) = [1, 2]

    sol = 1 + sin(cxy(1) * xy(1)) * sin(cxy(2) * xy(2)) * &
         exp(-sum(cxy**2) * diffusion_coeff * t)
  end function solution

  subroutine set_rhs(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(1:nc, 1:nc, i_rhs) = box%cc(1:nc, 1:nc, i_phi)
  end subroutine set_rhs

  subroutine prolong_to_new_children(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box(tree%boxes, id, i_phi, &
               a2_gc_interp_lim, a2_bc_dirichlet_zero)
       end do
    end do
  end subroutine prolong_to_new_children

  ! Compute L * phi, where L corresponds to (D * dt * nabla^2 - 1)
  subroutine box_op_diff(box, i_out, mg)
    type(box2_t), intent(inout) :: box    !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi
    real(dp)                    :: tmp

    nc    = box%n_cell
    tmp   = diffusion_coeff * dt / box%dr**2
    i_phi = mg%i_phi

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = tmp * ((4 + 1/tmp) * &
               box%cc(i, j, i_phi) &
               - box%cc(i-1, j, i_phi) &
               - box%cc(i+1, j, i_phi) &
               - box%cc(i, j-1, i_phi) &
               - box%cc(i, j+1, i_phi))
       end do
    end do
  end subroutine box_op_diff

  ! Locally solve L * phi = rhs, where L corresponds to (D * dt * nabla^2 - 1)
  subroutine box_gsrb_diff(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box            !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_rhs
    real(dp)                    :: tmp

    tmp   = diffusion_coeff * dt / box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = 1 / (4 + 1/tmp) * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) + &
               1/tmp * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine box_gsrb_diff

end program drift_diffusion_2d
