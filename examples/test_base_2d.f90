!> \example test_base_2d.f90
!> This example shows the basic functionality of m_afivo_2d.
program test_base
  use m_afivo_2d

  implicit none

  integer, parameter  :: dp           = kind(0.0d0)
  type(a2_t)          :: tree
  integer             :: i
  integer, parameter  :: n_boxes_base = 2
  integer             :: ix_list(2, n_boxes_base)
  integer             :: nb_list(4, n_boxes_base)
  integer, parameter  :: box_size     = 8
  integer, parameter  :: i_phi        = 1, i_mrtn = 2
  integer, parameter  :: n_var_cell   = 2
  integer, parameter  :: n_var_face   = 0
  integer, parameter  :: coarsen_to   = -1
  real(dp), parameter :: r_min(2)     = [0.0_dp, 0.0_dp]
  character(len=40)   :: var_names(2) = ["phi ", "mrtn"]
  real(dp)            :: dr
  character(len=40)   :: fname

  dr = 2 * acos(-1.0_dp) / box_size ! 2 * pi / box_size

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! Number of cells per coordinate in a box
       n_var_cell, &   ! Number of face-centered variables
       n_var_face, &   ! Number of cell-centered variables
       dr)             ! Distance between cells on base level

  ! Set up geometry
  ix_list(:, 1) = [1,1] ! One box at 1,1
  ix_list(:, 2) = [2,1] ! One box at 2,1

  ! Set neighbors for box one
  nb_list(a2_nb_lx, 1) = 2
  nb_list(a2_nb_hx, 1) = 2
  nb_list(a2_nb_ly, 1) = 1
  nb_list(a2_nb_hy, 1) = 1

  ! Set neighbors for box two
  nb_list(:, 2) = 2

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)

  ! Set variables on base
  call a2_loop_box(tree, set_init_cond)
  call a2_loop_boxes(tree, set_morton_variable)

  ! Fill ghost cells for phi
  call a2_gc_sides(tree, i_phi, a2_sides_prolong1, have_no_bc)

  do i = 1, 15
     print *, "i = ", i, "max_id", tree%max_id

     write(fname, "(A,I0,A)") "test_base_2d_", i, ".vtu"
     call a2_write_vtk(tree, trim(fname), var_names, i, i * 1.0_dp)

     call a2_adjust_refinement(tree, set_ref_flags)

     ! Tidy up will reorder all the boxes if reorder = .true.
     call a2_tidy_up(tree, max_frac_used=0.75_dp, goal_frac_used=0.5_dp, &
          n_clean_min=10000, reorder=.true.)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call a2_loop_boxes(tree, set_morton_variable)
  end do

  call a2_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    real(dp)                 :: rr

    call random_number(rr)
    if (rr < 0.2_dp .and. boxes(id)%lvl < 10) then
       ref_flags(id) = a5_do_ref
    else
       ref_flags(id) = a5_rm_ref
    end if
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_phi) = sin(0.5_dp * xy(1)) * cos(xy(2))
       end do
    end do
  end subroutine set_init_cond

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, i_phi, .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine set_morton_variable(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    boxes(id)%cc(:,:,i_mrtn) = id
  end subroutine set_morton_variable

  subroutine have_no_bc(boxes, id, i, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, iv
    stop "We have no boundary conditions in this example"
    boxes(id)%cc(1, i, iv) = 0    ! Prevent warning
  end subroutine have_no_bc

end program test_base
