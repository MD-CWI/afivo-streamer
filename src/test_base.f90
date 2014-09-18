program test_base
  use m_afivo

  implicit none

  type(a2_t) :: tree
  integer    :: i
  integer    :: ix(2)
  integer    :: ix_list(2, 1)
  integer    :: nb_list(4, 1)
  integer    :: box_size    = 16
  integer    :: n_boxes_max = 100
  real(dp)   :: dr(2)

  dr = 0.1_dp

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_cc=1, n_fx=0, n_fy=0)

  ! Create base level which is periodic in all directions
  ix               = (/1,1/)
  ix_list(:, 1)    = ix
  nb_list(:, 1) = 1

  call a2_set_base(tree, ix_list, nb_list, dr, (/0.0_dp, 0.0_dp/))
  call a2_loop_box(tree, set_init_cond)
  
  do i = 1, 1
     call a2_adjust_refinement(tree, ref_func)
  end do

  call a2_write_tree(tree, "test_base.silo", (/"my_var"/), &
       (/"my_unit"/), 0, 0.0_dp)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    if (all(box%ix == 1)) then
       ref_func = a5_do_ref
    else
       ref_func = a5_kp_ref
    end if
  end function ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, bs
    real(dp)                    :: x, y
    
    bs = size(box%cc, 1) - 2
    do j = 1, bs
       y = box%r_min(2) + (j-0.5_dp) * box%dr(2)
       do i = 1, bs
          x = box%r_min(1) + (i-0.5_dp) * box%dr(1)
          box%cc(i, j, 1) = sin(x * y)
       end do
    end do
  end subroutine set_init_cond
end program test_base
