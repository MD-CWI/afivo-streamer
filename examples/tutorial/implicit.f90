
program tutorial
  use m_a2_t
  use m_a2_core
  use m_a2_utils
  use m_a2_gc
  use m_a2_io
  use m_a2_mg

  implicit none

  integer, parameter  :: box_size     = 128
  integer, parameter  :: n_var_cell   = 3
  integer, parameter  :: n_var_face   = 0
  integer, parameter  :: i_phi        = 1
  integer, parameter  :: i_rhs        = 2
  integer, parameter  :: i_tmp        = 3

  real(dp), parameter :: domain_len   = 1.0_dp
  real(dp), parameter :: diff_coeff   = 1.0_dp
  real(dp), parameter :: dt_global    = 1.0e-2_dp

  type(a2_t)          :: tree
  type(mg2_t)         :: mg
  integer             :: output_cnt, n_steps
  real(dp)            :: dt_output, time, end_time
  character(len=100)  :: fname

  print *, "This is an Afivo tutorial for solving the heat equation."

  call initialize_tree(tree)
  call initialize_mg(mg)

  call a2_loop_box(tree, set_init_cond)
  call a2_gc_tree(tree, i_phi, a2_gc_interp_lim, a2_bc_dirichlet_zero)

  output_cnt = 0
  n_steps    = 0
  time       = 0.0_dp
  end_time   = 1.0_dp
  dt_output  = 0.02_dp

  do while (time < end_time)
     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "heat_implicit_", output_cnt
        call a2_write_vtk(tree, fname, ["phi", "rhs", "res"], output_cnt, time)
     end if

     call heat_backward_euler(tree, n_steps == 0)
     time = time + dt_global
     n_steps = n_steps + 1
  end do

contains

  subroutine initialize_tree(tree)
    type(a2_t), intent(inout) :: tree

    integer, parameter  :: n_boxes_base = 3
    real(dp)            :: dr
    integer             :: ix_list(2, n_boxes_base)
    integer             :: nb_list(4, n_boxes_base)

    dr = domain_len / box_size

    ! Initialize the tree
    call a2_init(tree, & ! Tree to initialize
         box_size, &     ! Number of cells per coordinate in a box
         n_var_cell, &   ! Number of face-centered variables
         n_var_face, &   ! Number of cell-centered variables
         dr, &           ! Distance between cells on base level
         coarsen_to=2)   ! Create extra coarse levels

    ! Set the indices for boxes 1 to 3
    ix_list(:, 1) = [1, 1]
    ix_list(:, 2) = [2, 1]
    ix_list(:, 3) = [2, 2]

    ! Negative values are used for boundary conditions
    nb_list(1:4, 1) = -1
    nb_list(1:4, 2) = -1
    nb_list(1:4, 3) = -1

    ! Connect the two boxes, nb_hx stands for "neighbor in high-x direction".
    ! Connections only have to be specified from one side.
    nb_list(a2_nb_hx, 1) = 2
    nb_list(a2_nb_hy, 2) = 3

    ! Construct the base (level 1) mesh
    call a2_set_base(tree, ix_list, nb_list)

    call a2_print_info(tree)
  end subroutine initialize_tree

  subroutine initialize_mg(mg)
    type(mg2_t), intent(inout) :: mg

    ! Set the multigrid options
    mg%i_phi        = i_phi
    mg%i_tmp        = i_tmp
    mg%i_rhs        = i_rhs
    mg%sides_bc     => a2_bc_dirichlet_zero
    mg%box_op       => lpl_heat
    mg%box_gsrb     => gsrb_heat

    call mg2_init_mg(mg)
  end subroutine initialize_mg

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy-0.5_dp) < 0.25) then
             box%cc(i, j, i_phi) = 1
          else
             box%cc(i, j, i_phi) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine heat_backward_euler(tree, first_time)
    type(a2_t), intent(inout) :: tree
    logical, intent(in)       :: first_time

    ! Set right-hand side
    call a2_loop_box(tree, set_rhs)

    ! Solve
    call mg2_fas_fmg(tree, mg, .true., first_time)
  end subroutine heat_backward_euler

  subroutine set_rhs(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    nc = box%n_cell
    box%cc(1:nc, 1:nc, i_rhs) = -box%cc(1:nc, 1:nc, i_phi) / dt_global
  end subroutine set_rhs

  !> Perform Laplacian operator on a box
  subroutine lpl_heat(box, i_out, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi
    real(dp)                    :: inv_dr_sq, fac

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    fac       = 4 + box%dr**2 / (dt_global * diff_coeff)

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_phi) + &
               box%cc(i+1, j, i_phi) + box%cc(i, j-1, i_phi) + &
               box%cc(i, j+1, i_phi) - fac * box%cc(i, j, i_phi))
       end do
    end do
  end subroutine lpl_heat

  !> Perform Gauss-Seidel relaxation on a box
  subroutine gsrb_heat(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_rhs
    real(dp)                    :: dx2, fac

    dx2   = box%dr**2
    fac   = 1 / (4 + dx2/(dt_global * diff_coeff))
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = fac * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine gsrb_heat

end program
