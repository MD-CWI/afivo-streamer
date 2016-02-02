!> \example test_mg_benchmark_2d.f90

! This program can be used to benchmark the multigrid routines. For simplicity,
! it does not compare results with known solution.
program test_mg2_2d
  use m_a2_t
  use m_a2_core
  use m_a2_mg
  use m_a2_utils
  use m_a2_gc
  use m_a2_io

  implicit none

  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 100
  integer, parameter :: n_var_cell = 3
  integer, parameter :: i_phi = 1
  integer, parameter :: i_tmp = 2
  integer, parameter :: i_rhs = 3

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  integer            :: n_cell, mesh_size, max_ref_lvl
  real(dp)           :: dr
  character(len=40)  :: arg_string
  type(mg2_t)        :: mg

  ! Get box size and mesh size from command line argument
  if (command_argument_count() /= 2) stop "Arguments should be: n_cell mesh_size"
  call get_command_argument(1, arg_string)
  read(arg_string, *) n_cell
  call get_command_argument(2, arg_string)
  read(arg_string, *) mesh_size

  ! Determine maximum refinement level
  max_ref_lvl = nint(log(mesh_size / real(n_cell, dp)) / log(2.0_dp)) + 1

  print *, "Box size: ", n_cell
  print *, "Mesh size:", 2**(max_ref_lvl-1) * n_cell

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of boxnn
  nb_list(:, 1) = -1            ! Dirichlet zero -> -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)

  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a2_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => a2_bc_dirichlet_zero ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  call mg2_init_mg(mg)

  ! Do the actual benchmarking
  do i = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg2_fas_fmg(tree, mg, .true., i>1)
  end do

  ! This writes a Silo output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  call a2_write_silo(tree, "poisson_benchmark_2d")

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a2_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(boxes, id, ref_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag

    ! Fully refine up to max_ref_lvl
    if (boxes(id)%lvl < max_ref_lvl) ref_flag = a5_do_ref
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(1:nc, 1:nc, i_rhs) = 1.0_dp
  end subroutine set_init_cond

end program test_mg2_2d
