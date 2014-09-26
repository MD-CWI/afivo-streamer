program test_diffusion
  use m_afivo

  implicit none

  type(a2_t)         :: tree
  integer            :: i, n, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 8
  integer            :: n_boxes_max = 100
  real(dp)           :: dr(2), dt(1)
  character(len=40)  :: fname

  dr    = 2 * acos(-1.0_dp) / box_size

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_cc=1, n_flux=1, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of box
  nb_list(:, id) = id            ! Box is periodic, so its own neighbor

  call a2_set_base(tree, ix_list, nb_list)
  do i = 1, 3
     call a2_adjust_refinement(tree, ref_func)
  end do
  call a2_loop_box(tree, set_init_cond)
  do i = 1, 3
     call a2_adjust_refinement(tree, ref_func)
  end do
  call a2_clear_tagbit(tree, a5_bit_new_children)

  call a2_loop_box(tree, set_init_cond)
  call a2_gc_sides(tree, a2_sides_from_parent)
  call a2_gc_corners(tree, a2_corners_from_parent)

  do i = 1, 20
     print *, "i = ", i, "n_boxes", tree%n_boxes

     write(fname, "(A,I0,A)") "test_diffusion_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), (/"my_var"/), i, i * 1.0_dp)

     do n = 1, 100
        ! Set timestep
        dt(1) = 0.5_dp * minval(a2_min_dr(tree))**2
        call a2_loop_box(tree, calculate_fluxes)
        call a2_consisent_fluxes(tree, [1])
        call a2_loop_box_arg(tree, update_solution, dt)
        call a2_gc_sides(tree, a2_sides_from_parent)
        call a2_gc_corners(tree, a2_corners_from_parent)
     end do

     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call a2_tidy_up(tree, 0.5_dp, 0.25_dp, 100*1000, .false.)
  end do

  call a2_destroy(tree)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    real(dp) :: max_dd(2)
    integer :: nc

    nc = box%cfg%n_cell
    max_dd(1) = maxval(abs(box%cc(0:nc-1, 1:nc, 1) - &
         2 * box%cc(1:nc, 1:nc, 1) + box%cc(2:nc+1, 1:nc, 1)))
    max_dd(2) = maxval(abs(box%cc(1:nc, 0:nc-1, 1) - &
         2 * box%cc(1:nc, 1:nc, 1) + box%cc(1:nc, 2:nc+1, 1)))
    max_dd = max_dd / a2_dr(box)**2

    if ((box%lvl < 3 .or. maxval(max_dd) > 0.6_dp) .and. box%lvl < 6) then
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

  subroutine calculate_fluxes(box)
    type(box2_t), intent(inout) :: box
    real(dp) :: inv_dr(2)
    integer :: nc

    nc = box%cfg%n_cell
    inv_dr = 1/a2_dr(box)

    box%fx(:,:,1) = box%cc(0:nc, 1:nc, 1) - box%cc(1:nc+1, 1:nc, 1)
    box%fx(:,:,1) = box%fx(:,:,1) * inv_dr(1)
    box%fy(:,:,1) = box%cc(1:nc, 0:nc, 1) - box%cc(1:nc, 1:nc+1, 1)
    box%fy(:,:,1) = box%fy(:,:,1) * inv_dr(2)
  end subroutine calculate_fluxes

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr(2)
    integer                     :: nc

    nc = box%cfg%n_cell
    inv_dr = 1/a2_dr(box)

    nc = box%cfg%n_cell
    box%cc(1:nc, 1:nc, 1) = box%cc(1:nc, 1:nc, 1) + dt(1) * ( &
         (box%fx(1:nc, :, 1) - box%fx(2:nc+1, :, 1)) * inv_dr(1) + &
         (box%fy(:, 1:nc, 1) - box%fy(:, 2:nc+1, 1)) * inv_dr(2))
  end subroutine update_solution

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, [1], .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

end program test_diffusion
