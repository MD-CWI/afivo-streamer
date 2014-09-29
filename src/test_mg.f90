program test_mg
  use m_afivo
  use m_mg2d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_res = 4
  type(a2_t)         :: tree
  integer            :: i, id, n_changes
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 2
  integer            :: n_boxes_max = 10*1000
  real(dp)           :: dr(2)
  character(len=40)  :: fname, var_names(4)

  type(mg2_t) :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_res) = "res"
  call mg2d_set(mg, i_phi, i_tmp, i_rhs, i_res, 2, 2, &
       sides_bc, a2_corners_extrap)

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_var_cell=4, n_var_face=0, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do i = 1, 20
     call a2_adjust_refinement(tree, ref_func_init, n_changes)
     if (n_changes == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  ! Restrict the rhs from children recursively
  call mg2d_restrict_rhs(tree, mg)

  do i = 1, 10
     write(fname, "(A,I0,A)") "test_mg_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), var_names, i, 0.0_dp)

     ! call mg2d_fas_vcycle(tree, mg, tree%n_lvls)
     call mg2d_fas_fmg(tree, mg)
  end do

  print *, "n_boxes", tree%n_boxes
  print *, "n_cells", tree%n_boxes * tree%cfg%n_cell**2

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    if (box%lvl < 4 .or. &
         (box%lvl < 7 .and. (norm2(a2_r_center(box)-2) < 1.25_dp))) then
       ref_func_init = a5_do_ref
    else
       ref_func_init = a5_rm_ref
    end if
  end function ref_func_init

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%cfg%n_cell
    box%cc(:, :, i_phi) = 0

    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 1) then
             box%cc(i, j, i_rhs) = 3
          else
             box%cc(i, j, i_rhs) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine sides_bc(boxes, id, nb, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, ivs(:)
    integer                     :: nc
    real(dp) :: dr(2)

    nc = boxes(id)%cfg%n_cell

    if (boxes(id)%neighbors(nb) == -1) then
       select case (nb)
       case (a2_nb_lx)
          ! Dirichlet zero
          boxes(id)%cc(0, 1:nc, ivs) = -boxes(id)%cc(1, 1:nc, ivs)
       case (a2_nb_hx)
          ! Derivative of 2
          dr = a2_dr(boxes(id))
          boxes(id)%cc(nc+1, 1:nc, ivs) = 2 * dr(1) + &
               boxes(id)%cc(nc, 1:nc, ivs)
       case (a2_nb_ly)
          ! Neumann zero
          boxes(id)%cc(1:nc, 0, ivs) = boxes(id)%cc(1:nc, 1, ivs)
       case (a2_nb_hy)
          ! Neumann zero
          boxes(id)%cc(1:nc, nc+1, ivs) = boxes(id)%cc(1:nc, nc, ivs)
       end select
    end if
  end subroutine sides_bc

end program test_mg