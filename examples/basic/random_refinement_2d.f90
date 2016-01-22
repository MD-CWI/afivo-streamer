!> \example random_refinement_2d.f90

! This example shows how to create an AMR tree, perform random refinement, write
! output files, and how to fill ghost cells.
program random_refinement_2d
  use m_a2_t
  use m_a2_core
  use m_a2_utils
  use m_a2_io
  use m_a2_gc

  implicit none

  type(a2_t)           :: tree
  integer              :: i
  integer, parameter   :: n_boxes_base = 1
  integer              :: ix_list(2, n_boxes_base)
  integer              :: nb_list(4, n_boxes_base)
  integer, parameter   :: n_cell       = 8
  integer, parameter   :: i_phi        = 1
  integer, parameter   :: n_var_cell   = 2
  integer, parameter   :: n_var_face   = 0
  type(ref_info_t)     :: ref_info
  real(dp)             :: dr
  character(len=40)    :: fname

  ! The cell spacing at the coarsest grid level
  dr = 2 * acos(-1.0_dp) / n_cell ! 2 * pi / n_cell

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       1, &            ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       cc_names = ["phi"])      ! Optional: names of cell-centered variables

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1] ! One box at index 1,1

  ! Set neighbors for box one, here nb means neighbor (direction) and l/h stands
  ! for low/high
  nb_list(a2_nb_lx, 1) = 1      ! lower-x neighbor of box 1 is box 1
  nb_list(a2_nb_hx, 1) = 1      ! higher-x neighbor of box 1 is box 1
  nb_list(a2_nb_ly, 1) = 1      ! etc. for y-direction
  nb_list(a2_nb_hy, 1) = 1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)

  ! Set variables on base by using the helper functions a2_loop_box(tree, sub)
  ! and a2_loop_boxes(tree, sub). These functions call the subroutine sub for
  ! each box in the tree, with a slightly different syntax.
  call a2_loop_box(tree, set_init_cond)

  ! Fill ghost cells for phi. The third argument is a subroutine that fills
  ! ghost cells near refinement boundaries, and the fourth argument fill ghost
  ! cells near physical boundaries.
  call a2_gc_tree(tree, i_phi, a2_gc_interp, a2_bc_dirichlet_zero)

  do i = 1, 10
     ! This writes a VTK output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     write(fname, "(A,I0)") "random_refinement_2d_", i
     call a2_write_vtk(tree, trim(fname), dir="output")

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     print *, "# new     boxes", ref_info%n_add
     print *, "# removed boxes", ref_info%n_rm

     ! Newly added boxes should be initialized, which is done in this routine.
     call prolong_to_new_children(tree, ref_info)
  end do

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a2_destroy(tree)

contains

  ! This routine sets refinement flags
  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:) ! A list of all boxes in the tree
    integer, intent(in)      :: id       ! The index of the current box
    integer, intent(inout)   :: ref_flags(:) ! A list with the refinement flags
                                             ! for the boxes
    real(dp)                 :: rr

    ! Draw a [0, 1) random number
    call random_number(rr)

    if (rr < 0.25_dp .and. boxes(id)%lvl < 8) then
       ref_flags(id) = a5_do_ref ! Add refinement
    else if (rr > 0.75_dp) then
       ref_flags(id) = a5_rm_ref ! Ask to remove this box, which will not always
                                 ! happen (see documentation)
    else
       ref_flags(id) = a5_kp_ref ! Keep the box as-is (which is the default
                                 ! action if you don't specify anything)
    end if
  end subroutine set_ref_flags

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          ! Get the coordinate of the cell center at i,j
          xy = a2_r_cc(box, [i,j])

          ! Set the values at each cell according to some function
          box%cc(i, j, i_phi) = sin(0.5_dp * xy(1)) * cos(xy(2))
       end do
    end do
  end subroutine set_init_cond

  ! Set values on newly added boxes
  subroutine prolong_to_new_children(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       ! For each newly added box ...
       do i = 1, size(ref_info%lvls(lvl)%add)
          ! Use prolongation to set the values of variable i_phi
          id = ref_info%lvls(lvl)%add(i)
          call a2_prolong2_to(tree%boxes, id, i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          ! After values have been set on this level, fill ghost cells
          call a2_gc_box(tree%boxes, id, i_phi, &
               a2_gc_interp, a2_bc_dirichlet_zero)
       end do
    end do
  end subroutine prolong_to_new_children

end program random_refinement_2d
