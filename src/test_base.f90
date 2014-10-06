program test_base
  use m_afivo

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  type(a2_t)         :: tree
  integer            :: i, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 8
  integer, parameter :: i_phi = 1, i_mrtn = 2
  character(len=40)  :: var_names(2) = ["phi ", "mrtn"]
  integer            :: n_boxes_max = 100
  integer            :: n_lvls_max = 20
  real(dp)           :: dr
  character(len=40)  :: fname

  dr = 2 * acos(-1.0_dp) / box_size

  ! Initialize tree
  call a2_init(tree, n_lvls_max, n_boxes_max, box_size, n_var_cell=2, &
       n_var_face=0, dr = dr, r_min = [0.0_dp, 0.0_dp])

  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = id    ! Box is periodic, so its own neighbor

  call a2_set_base(tree, ix_list, nb_list)

  ! Set variables on base
  call a2_loop_box(tree, set_init_cond)
  call a2_loop_boxes(tree, set_morton_variable)

  ! Fill ghost cells for phi
  call a2_gc_sides(tree, [i_phi], a2_sides_prolong1, have_no_bc)
  call a2_gc_corners(tree, [i_phi], a2_corners_prolong1, have_no_bc)

  do i = 1, 20
     print *, "i = ", i, "max_id", tree%max_id

     write(fname, "(A,I0,A)") "test_base_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), var_names, i, i * 1.0_dp)

     call a2_adjust_refinement(tree, ref_func)
     call a2_tidy_up(tree, max_frac_used=0.75_dp, goal_frac_used=0.5_dp, &
          n_clean_min=10000, just_reorder=.true.)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call a2_loop_boxes(tree, set_morton_variable)
  end do

  call a2_destroy(tree)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    real(dp)                 :: rr

    call random_number(rr)
    if (rr < 0.2_dp) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
  end function ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_phi) = sin(xy(1)) * cos(xy(2))
       end do
    end do
  end subroutine set_init_cond

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, [i_phi], .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine set_morton_variable(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    boxes(id)%cc(:,:,i_mrtn) = id
  end subroutine set_morton_variable

  subroutine have_no_bc(boxes, id, i, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, ivs(:)
    stop "We have no boundary conditions in this example"
  end subroutine have_no_bc

end program test_base
