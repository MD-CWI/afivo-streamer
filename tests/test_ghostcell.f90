! Test the refinement procedure
program test_init
  use m_a2_types
  use m_a2_core
  use m_a2_ghostcell
  use m_a2_utils

  use m_a3_types
  use m_a3_core
  use m_a3_ghostcell
  use m_a3_utils

  implicit none

  type(a2_t)       :: tree_2d
  type(a3_t)       :: tree_3d
  type(ref_info_t) :: ref_info
  integer          :: i, n_lvl
  integer          :: ixs_2d(2, 1), nbs_2d(a2_num_neighbors, 1)
  integer          :: ixs_3d(3, 1), nbs_3d(a3_num_neighbors, 1)

  call a2_init(tree_2d, n_cell=8, n_var_cell=1, n_var_face=0, dr = 1.0_dp)
  call a3_init(tree_3d, n_cell=8, n_var_cell=1, n_var_face=0, dr = 1.0_dp)

  ixs_2d = 1                    ! Box at 1,1
  nbs_2d = -1                   ! Boundary condition indicated by -1
  call a2_set_base(tree_2d, ixs_2d, nbs_2d)

  ixs_3d = 1                    ! Box at 1,1,1
  nbs_3d = -1                   ! Boundary condition indicated by -1
  call a3_set_base(tree_3d, ixs_3d, nbs_3d)

  n_lvl = 4

  do i = 1, n_lvl
     call a2_adjust_refinement(tree_2d, refinement_2d, ref_info)
     call a3_adjust_refinement(tree_3d, refinement_3d, ref_info)
  end do

  call a2_loop_box(tree_2d, init_2d)
  call a3_loop_box(tree_3d, init_3d)

  ! Should set all ghost cells to zero
  call a2_gc_tree(tree_2d, 1, a2_gc_interp, a2_bc_neumann_zero)
  call a3_gc_tree(tree_3d, 1, a3_gc_interp, a3_bc_neumann_zero)

  call a2_loop_box(tree_2d, check_ghostcell_2d)
  call a3_loop_box(tree_3d, check_ghostcell_3d)

  print *, "Success"

contains

  subroutine refinement_2d(boxes, id, ref_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag

    if (all(boxes(id)%ix == 1)) ref_flag = af_do_ref
  end subroutine refinement_2d

  subroutine refinement_3d(boxes, id, ref_flag)
    type(box3_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag

    if (all(boxes(id)%ix == 1)) ref_flag = af_do_ref
  end subroutine refinement_3d

  subroutine init_2d(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc

    nc                    = box%n_cell
    box%cc(1:nc, 1:nc, 1) = 0
  end subroutine init_2d

  subroutine init_3d(box)
    type(box3_t), intent(inout) :: box
    integer                     :: nc

    nc                    = box%n_cell
    box%cc(1:nc, 1:nc, 1:nc, 1) = 0
  end subroutine init_3d

  subroutine check_ghostcell_2d(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: tmp

    nc  = box%n_cell
    tmp = sum(box%cc(1:nc, 0, 1)) + &
         sum(box%cc(1:nc, nc+1, 1)) + &
         sum(box%cc(0, 1:nc, 1)) + &
         sum(box%cc(nc+1, 1:nc, 1))

    if (tmp > 0 .or. tmp < 0) stop "Wrong ghostcell value"
  end subroutine check_ghostcell_2d

  subroutine check_ghostcell_3d(box)
    type(box3_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: tmp

    nc  = box%n_cell
    tmp = sum(box%cc(1:nc, 1:nc, 0, 1)) + &
         sum(box%cc(1:nc, 1:nc, nc+1, 1)) + &
         sum(box%cc(1:nc, 0, 1:nc, 1)) + &
         sum(box%cc(1:nc, nc+1, 1:nc, 1)) + &
         sum(box%cc(0, 1:nc, 1:nc, 1)) + &
         sum(box%cc(nc+1, 1:nc, 1:nc, 1))

    if (tmp > 0 .or. tmp < 0) stop "Wrong ghostcell value"
  end subroutine check_ghostcell_3d

end program
