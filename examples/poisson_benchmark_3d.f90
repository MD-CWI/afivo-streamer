!> \example test_mg_benchmark_3d.f90

! This program can be used to benchmark the multigrid routines. For simplicity,
! it does not compare results with known solution.
program poisson_benchmark_3d
  use m_a3_t
  use m_a3_core
  use m_a3_mg
  use m_a3_utils
  use m_a3_gc
  use m_a3_io

  implicit none

  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 10
  integer, parameter :: n_var_cell = 3
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_tmp = 3

  type(a3_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: box_size, mg_iter
  integer            :: ix_list(3, n_boxes_base)
  integer            :: nb_list(6, n_boxes_base)
  integer            :: mesh_size, max_ref_lvl
  real(dp)           :: dr
  character(len=40)  :: fname,filename
  type(mg3_t)        :: mg
  integer            :: count_rate,t_start, t_end

  write(*,'(A)') 'program poisson_benchmark_3d'

  call parallel_threads()

  filename='input'
  ! Get box size and mesh size from file input_benchmark
  open(unit=5,file=trim(filename),status='UNKNOWN',position='REWIND')
  read(5, *) box_size
  read(5, *) mesh_size

  ! Determine maximum refinement level
  max_ref_lvl = nint(log(mesh_size / real(box_size, dp)) / log(2.0_dp)) + 1

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  write(*,'(A,i4)') 'Box size             :', box_size
  write(*,'(A,i4)') 'Mesh size            :', 2**(max_ref_lvl-1) * box_size
  write(*,'(A,i4)') 'Max refinement level :', max_ref_lvl

  ! Initialize tree
  call a3_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1,1]       ! Set index of boxnn

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1            ! Dirichlet zero -> -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a3_set_base(tree, ix_list, nb_list)
  call a3_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call a3_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a3_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)
  write(*,'(A,f8.2,1x,A,1/)') 'Time making amr grid = ', &
          (t_end-t_start) / real(count_rate,dp), &
          ' seconds'

  call a3_print_info(tree)
  dr = a3_min_dr(tree)
  write(*,'(A,2x,Es11.4)') ' dr of finest level:',dr

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => a3_bc_dirichlet_zero ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call
  ! of mg3_fas_fmg
  call mg3_init_mg(mg)

  ! Do the actual benchmarking
  call system_clock(t_start, count_rate)
  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg3_fas_fmg(tree, mg, .true., mg_iter>1)

     ! This writes a VTK output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     write(fname, "(A,I0)") "poisson_benchmark_3d_", mg_iter
     call a3_write_vtk(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, '(A,i3,1x,A,f8.2,1x,A,/)') &
           ' Wall-clock time after ',n_iterations, &
           ' multigrid iterations: ', (t_end-t_start) / real(count_rate, dp), &
           ' seconds'
   
  ! This writes a VTK output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  call a3_write_vtk(tree, "poisson_benchmark_3d", dir="output")

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a3_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(boxes, id, ref_flag)
    type(box3_t), intent(in) :: boxes(:) ! A list of all boxes in the tree
    integer, intent(in)      :: id       ! The index of the current box
    integer, intent(inout)   :: ref_flag

    ! Fully refine up to max_ref_lvl
    if (boxes(id)%lvl < max_ref_lvl) ref_flag = a5_do_ref
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box3_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(1:nc, 1:nc, 1:nc, i_rhs) = 1.0_dp
  end subroutine set_init_cond

end program poisson_benchmark_3d
