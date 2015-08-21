!> \example test_mg_2d.f90
!> Example showing how to use m_mg_2d
program test_mg
  use m_afivo_2d
  use m_mg_2d

  implicit none

  integer, parameter :: dp           = kind(0.0d0)
  integer, parameter :: box_size     = 4
  integer, parameter :: n_boxes_base = 3
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, id
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(3)
  type(mg2_t)        :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=3, n_var_face=0, dr = dr, &
       n_boxes = 10*1000)

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1
  nb_list(a2_nb_hx, id) = 2

  id = 2
  ix_list(:, id) = [2,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1
  nb_list(a2_nb_lx, id) = 1

  id = 3
  ix_list(:, id) = [2,2]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1
  nb_list(a2_nb_ly, id) = 2

  call a2_set_base(tree, ix_list, nb_list)

  do i = 1, 20
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  ! Set the multigrid options
  mg%i_phi        = i_phi
  mg%i_tmp        = i_tmp
  mg%i_rhs        = i_rhs
  mg%sides_bc     => sides_bc

  call mg2_init_mg(mg)

  do i = 1, 100
     ! call mg2_fas_vcycle(tree, mg, tree%max_lvl)
     call mg2_fas_fmg(tree, mg, .false., i == 1)
  end do

  write(fname, "(A,I0)") "test_mg_2d_", 1
  call a2_write_silo(tree, trim(fname), var_names, 1, 0.0_dp)

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)

    if (boxes(id)%lvl < 6 .or. &
         (boxes(id)%lvl < 8 .and. &
         (norm2(a2_r_center(boxes(id))-2) < 0.75_dp))) then
       ref_flags(id) = a5_do_ref
    else
       ref_flags(id) = a5_rm_ref
    end if
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    box%cc(:, :, i_phi) = 0

    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_rhs) = exp(-sum((xy - 2)**2))
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
          ! Dirichlet zero
          boxes(id)%cc(0, 1:nc, iv) = -boxes(id)%cc(1, 1:nc, iv)
       case (a2_nb_hx)
          ! Dirichlet zero
          boxes(id)%cc(nc+1, 1:nc, iv) = -boxes(id)%cc(nc, 1:nc, iv)
       case (a2_nb_ly)
          ! Dirichlet zero
          boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
       case (a2_nb_hy)
          ! Dirichlet zero
          boxes(id)%cc(1:nc, nc+1, iv) = -boxes(id)%cc(1:nc, nc, iv)
       end select
    end if
  end subroutine sides_bc

end program test_mg
