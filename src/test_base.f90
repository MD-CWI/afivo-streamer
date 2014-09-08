program test_base
  use m_afivo
  
  implicit none

  type(a2_t) :: tree

  ! Initialize tree
  call a2_init(tree, dr, n_boxes_max, box_size, n_cc=1)
  
  ! Create base level
  ix = (/1,1/)
  call a2_add_base_box(tree, ix, lnb=ix, rnb=ix, bnb=ix, tnb=ix)
  call a2_adjust_ref(tree, ref_func)
  call a2_write_tree(tree, 
end program test_base
