program test_base
  use m_afivo

  implicit none

  type(a2_t) :: tree
  integer    :: i
  integer    :: ix(2)
  integer    :: ix_list(2, 1)
  integer    :: nb_list(4, 1)
  integer, parameter :: box_size    = 8
  integer    :: n_boxes_max = 1000*1000
  real(dp)   :: dr(2)
  character(len=40) :: fname

  dr = 2 * acos(-1.0_dp) / box_size

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_cc=1, n_fx=0, n_fy=0, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  ! Create base level which is periodic in all directions
  ix            = (/1,1/)
  ix_list(:, 1) = ix
  nb_list(:, 1) = 1

  call a2_set_base(tree, ix_list, nb_list)
  call a2_loop_box(tree, set_init_cond)

  do i = 1, 25
     print *, "Iteration", i
     write(fname, "(A,I0,A)") "test_base_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), (/"my_var"/), i, i * 1.0_dp)

     call a2_fill_gc_sides(tree, gc_side_subr)
     call a2_fill_gc_corners(tree, gc_corner_subr)

     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)

     call a5_tidy_storage(tree, n_boxes_max)
  end do

  print *, "n_grids", tree%n_boxes
  call a2_destroy(tree)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    real(dp) :: rr
    call random_number(rr)

    if (rr < 0.2_dp) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
    ! if (all(box%ix == 1)) then
    !    ref_func = a5_do_ref
    ! else
    !    ref_func = a5_kp_ref
    ! end if
  end function ref_func

  subroutine gc_side_subr(boxes, id, nb)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb
    integer                     :: n, nc

    nc    = boxes(id)%cfg%n_cell
    select case (nb)
    case (nb_lx)
       call a2_prolong1_to(boxes, id, [0], [(n, n=1,nc)], [1])
    case (nb_hx)
       call a2_prolong1_to(boxes, id, [nc+1], [(n, n=1,nc)], [1])
    case (nb_ly)
       call a2_prolong1_to(boxes, id, [(n, n=1,nc)], [0], [1])
    case (nb_hy)
       call a2_prolong1_to(boxes, id, [(n, n=1,nc)], [nc+1], [1])
    end select
  end subroutine gc_side_subr

  subroutine gc_corner_subr(boxes, id, cn)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, cn
    integer                     :: nc

    nc    = boxes(id)%cfg%n_cell
    select case (cn)
    case (cn_lxly)
       call a2_prolong1_to(boxes, id, [0], [0], [1])
    case (cn_hxly)
       call a2_prolong1_to(boxes, id, [nc+1], [0], [1])
    case (cn_lxhy)
       call a2_prolong1_to(boxes, id, [0], [nc+1], [1])
    case (cn_hxhy)
       call a2_prolong1_to(boxes, id, [nc+1], [nc+1], [1])
    end select
  end subroutine gc_corner_subr

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, bs
    real(dp)                    :: xy(2)

    bs = size(box%cc, 1) - 2
    do j = 1, bs
       do i = 1, bs
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, 1) = sin(xy(1)) * cos(xy(2))
       end do
    end do
  end subroutine set_init_cond

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, [1])
    end if
    boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
  end subroutine prolong_to_new_children
end program test_base
