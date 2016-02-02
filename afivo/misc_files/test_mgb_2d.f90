!> \example test_mgb_2d.f90
!> Example showing how to use a level set function for boundary conditions
program test_mgb
  use m_a2_t
  use m_a2_core
  use m_a2_mg
  use m_a2_utils
  use m_a2_io

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_lsf = 4
  integer, parameter :: i_bval = 5
  real(dp), parameter :: phi0 = 1e3_dp

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, id
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(5)
  type(mg2_t)       :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_lsf) = "lsf"
  var_names(i_bval) = "bval"

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=5, n_var_face=0, &
       dr = dr, coarsen_to = 8, n_boxes = 10*1000)

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of box
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do i = 1, 20
     call a2_adjust_refinement(tree, ref_func, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  ! Set the multigrid options
  mg%i_phi       = i_phi
  mg%i_tmp       = i_tmp
  mg%i_rhs       = i_rhs
  mg%i_lsf       = i_lsf
  mg%i_bval      = i_bval
  mg%sides_bc    => sides_bc
  mg%box_op      => mg2_auto_op
  mg%box_corr    => mg2_auto_corr
  mg%box_gsrb    => mg2_auto_gsrb

  call mg2_init_mg(mg)

  do i = 1, 20
     ! call mg2_fas_vcycle(tree, mg, tree%highest_lvl, .true.)
     call mg2_fas_fmg(tree, mg, .true., i == 1)
     write(fname, "(A,I0)") "test_mgb_2d_", i
     call a2_write_vtk(tree, trim(fname), var_names, i, 0.0_dp)
  end do

  print *, "highest_id", tree%highest_id
  print *, "n_cells", tree%highest_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  subroutine ref_func(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)

    if (boxes(id)%lvl < 6) &
         ref_flags(id) = a5_do_ref
  end subroutine ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    box%cc(:, :, i_phi) = 0

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j]) - [2, 2]
          box%cc(i, j, i_rhs) = 0
          ! box%cc(i, j, i_lsf) = (sum(xy**2)-1)**3 - xy(1)**2 * xy(2)**3 ! Heart
          ! box%cc(i, j, i_lsf) = norm2(xy)-1.0_dp ! Circle
          box%cc(i, j, i_lsf) = sqrt(100 * xy(1)**2 + xy(2)**2) - 1.8_dp
          box%cc(i, j, i_bval) = phi0
       end do
    end do
  end subroutine set_init_cond

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_neighb_lowx)
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_neighb_highx)
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_neighb_lowy)
       boxes(id)%cc(1:nc, 0, iv) = - boxes(id)%cc(1:nc, 1, iv)
    case (a2_neighb_highy)
       boxes(id)%cc(1:nc, nc+1, iv) = - boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc

end program
