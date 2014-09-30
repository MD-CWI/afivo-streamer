program test_base
  use m_afivo

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  type(a2_t)         :: tree
  integer            :: i, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 8
  integer            :: n_boxes_max = 100
  real(dp)           :: dr
  character(len=40)  :: fname

  dr = 2 * acos(-1.0_dp) / box_size

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_var_cell=1, n_var_face=0, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of box
  nb_list(:, id) = id            ! Box is periodic, so its own neighbor

  call a2_set_base(tree, ix_list, nb_list)
  call a2_loop_box(tree, set_init_cond)
  call a2_gc_sides(tree, [1], a2_sides_prolong1, have_no_bc)
  call a2_gc_corners(tree, [1], a2_corners_prolong1, have_no_bc)

  do i = 1, 20
     print *, "i = ", i, "n_boxes", tree%n_boxes

     write(fname, "(A,I0,A)") "test_base_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), (/"my_var"/), i, i * 1.0_dp)

     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call a2_tidy_up(tree, 0.5_dp, 0.25_dp, 100*1000, .false.)
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

    nc = box%cfg%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, 1) = sin(xy(1)) * cos(xy(2))
       end do
    end do
  end subroutine set_init_cond

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, [1], .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine have_no_bc(boxes, id, i, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, ivs(:)
    stop "We have no boundary conditions in this example"
  end subroutine have_no_bc

end program test_base
