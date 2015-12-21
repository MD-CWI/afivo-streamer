!> \example test_base_3d.f90
!> Shows the basic functionality of m_a3_t
program test_base
  use m_a3_t
  use m_a3_core
  use m_a3_utils
  use m_a3_gc
  use m_a3_io
  use m_a3_prolong

  implicit none

  type(a3_t)         :: tree
  integer            :: i, id
  integer            :: ix_list(3, 1)
  integer            :: nb_list(6, 1)
  integer, parameter :: box_size    = 8
  integer, parameter :: i_phi = 1, i_mrtn = 2
  type(ref_info_t)   :: ref_info
  character(len=40)  :: var_names(2) = ["phi ", "mrtn"]
  real(dp)           :: dr
  character(len=40)  :: fname

  dr = 2 * acos(-1.0_dp) / box_size

  ! Initialize tree
  call a3_init(tree, box_size, n_var_cell=2, n_var_face=0, dr = dr)

  id             = 1
  ix_list(:, id) = [1,1,1]      ! Set spatial index of box
  nb_list(:, id) = id           ! Box is periodic, so its own neighbor

  call a3_set_base(tree, ix_list, nb_list)

  ! Set variables on base
  call a3_loop_box(tree, set_init_cond)
  call a3_loop_boxes(tree, set_morton_variable)

  ! Fill ghost cells for phi
  call a3_gc_tree(tree, i_phi, a3_gc_interp, have_no_bc)

  do i = 1, 13
     write(fname, "(A,I0)") "test_base_3d_", i
     call a3_write_vtk(tree, trim(fname), var_names, i, i * 1.0_dp)

     call a3_adjust_refinement(tree, set_ref_flags, ref_info)
     call prolong_to_new_children(tree, ref_info)

     call a3_tidy_up(tree, max_frac_used=0.75_dp, goal_frac_used=0.5_dp, &
          n_clean_min=10000, reorder=.true.)

     call a3_loop_boxes(tree, set_morton_variable)
  end do

  call a3_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box3_t), intent(in) :: boxes(:)
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
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             xyz = a3_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_phi) = sin(xyz(1)) * sin(xyz(2)) * sin(xyz(3))
          end do
       end do
    end do
  end subroutine set_init_cond

  subroutine prolong_to_new_children(tree, ref_info)
    type(a3_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a3_prolong1_to(tree%boxes, id, i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a3_gc_box(tree%boxes, id, i_phi, &
               a3_gc_interp, have_no_bc)
       end do
    end do
  end subroutine prolong_to_new_children

  subroutine set_morton_variable(boxes, id)
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    boxes(id)%cc(:,:,:,i_mrtn) = id
  end subroutine set_morton_variable

  subroutine have_no_bc(boxes, id, i, iv)
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, iv
    stop "We have no boundary conditions in this example"
    boxes(id)%cc(1, 1, i, iv) = 0    ! Prevent warning unused
  end subroutine have_no_bc

end program test_base
