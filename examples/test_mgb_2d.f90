!> \example test_mgb_2d.f90
!> Example showing how to use m_mgb_2d
program test_mgb
  use m_afivo_2d
  use m_mgb_2d

  implicit none

  integer, parameter :: dp           = kind(0.0d0)
  integer, parameter :: box_size     = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_res = 4
  integer, parameter :: i_lsf = 5

  type(a2_t)         :: tree
  integer            :: i, id, n_changes
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(5)
  type(mg2_t)       :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_res) = "res"
  var_names(i_lsf) = "lsf"

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=5, n_var_face=0, dr = dr)

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of box
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do i = 1, 20
     call a2_adjust_refinement(tree, ref_func_init, n_changes)
     if (n_changes == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  ! Set the multigrid options
  call mg2_set(mg, i_phi, i_tmp, i_rhs, i_res, i_lsf, 2, 2, 3, &
       sides_bc, mg2_lpl_box, mg2_gsrb_lpl_box)

  ! Restrict from children recursively
  call a2_restrict_tree(tree, i_rhs)
  call a2_restrict_tree(tree, i_phi)

  do i = 1, 20
     ! call mg2_fas_vcycle(tree, mg, tree%max_lvl)
     call mg2_fas_fmg(tree, mg)
     write(fname, "(A,I0,A)") "test_mgb_2d_", i, ".vtu"
     call a2_write_vtk(tree, trim(fname), var_names, i, 0.0_dp)
  end do

  ! write(fname, "(A,I0,A)") "test_mg_", 1, ".vtu"
  ! call a2_write_vtk(tree, trim(fname), var_names, 1, 0.0_dp)

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    if (box%lvl < 6) then
       ref_func_init = a5_do_ref
    else
       ref_func_init = a5_rm_ref
    end if
  end function ref_func_init

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    box%cc(:, :, i_phi) = 0

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j]) - 2
          box%cc(i, j, i_rhs) = 0 * exp(-sum((xy)**2))
          box%cc(i, j, i_lsf) = (sum(xy**2)-1)**3 - xy(1)**2 * xy(2)**3
       end do
    end do
  end subroutine set_init_cond

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    if (boxes(id)%neighbors(nb) == -1) then
       select case (nb)
       case (a2_nb_lx)
          boxes(id)%cc(0, 1:nc, iv) = 2-boxes(id)%cc(1, 1:nc, iv)
       case (a2_nb_hx)
          boxes(id)%cc(nc+1, 1:nc, iv) = 2-boxes(id)%cc(nc, 1:nc, iv)
       case (a2_nb_ly)
          boxes(id)%cc(1:nc, 0, iv) = 2-boxes(id)%cc(1:nc, 1, iv)
       case (a2_nb_hy)
          boxes(id)%cc(1:nc, nc+1, iv) = 2-boxes(id)%cc(1:nc, nc, iv)
       end select
    end if
  end subroutine sides_bc

end program test_mgb
