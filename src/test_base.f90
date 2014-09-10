program test_base
  use m_afivo
  
  implicit none

  type(a2_t) :: tree
  integer :: ix_list(2, 1)
  integer :: nb_list(2, 2, 1)

  ! Initialize tree
  call a2_init(tree, dr, n_boxes_max, box_size, n_cc=1)
  
  ! Create base level which is periodic in all directions
  ix               = (/1,1/)
  ix_list(:, 1)    = ix
  nb_list(:, :, 1) = 1
  
  call a2_set_base(tree, ix_list, nb_list)
  call a2_callback(tree, set_init_cond)
  call a2_adjust_ref(tree, ref_func)
  call a2_write_tree(tree, "test_base.silo", (/"my_var"/), &
       (/"my_unit"/), 0, 0.0_dp)
end program test_base
