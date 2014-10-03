! NOTE: this is a test thingy and it is not working
program test_mg
  use m_afivo
  use m_mg2d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_res = 4
  integer, parameter :: i_eta = 5, i_pot = 6
  type(a2_t)         :: tree
  integer            :: i, id, n_changes
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 16
  integer            :: n_boxes_max = 10*1000
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(6)

  type(mg2_t)      :: mg
  type(mg2_subt_t) :: subtree

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_res) = "res"
  var_names(i_eta) = "eta"
  var_names(i_pot) = "pot"

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, 20, n_boxes_max, box_size, n_var_cell=6, n_var_face=0, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do i = 1, 0
     call a2_adjust_refinement(tree, ref_func_init, n_changes)
     if (n_changes == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  call mg2d_set(mg, i_phi, i_tmp, i_rhs, i_res, 2, 2, 2, &
       sides_bc, a2_corners_extrap, elec_lpl, elec_gsrb)

  call mg2d_create_subtree(tree, subtree)

  ! Restrict from children recursively
  call mg2d_restrict_trees(tree, subtree, [i_rhs, i_phi, i_eta, i_pot], mg, .true.)

  do i = 1, 10
     ! call mg2d_fas_vcycle(tree, subtree, mg, .true., tree%n_lvls)
     call mg2d_fas_fmg(tree, subtree, mg, .false.)
     write(fname, "(A,I0,A)") "test_mg_elec_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), var_names, i, 0.0_dp)
  end do

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    if (box%lvl < 4 .or. &
         (box%lvl < 6 .and. (norm2(a2_r_center(box)-2) < 1.0_dp))) then
       ref_func_init = a5_do_ref
    else
       ref_func_init = a5_rm_ref
    end if
  end function ref_func_init

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), eta, pot, inv_dr_sq

    nc = box%n_cell
    inv_dr_sq = 1/box%dr**2
    box%cc(:, :, i_phi) = 0

    do j = 1, nc
       do i = 1, nc

          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 0.75_dp) then
             box%cc(i, j, i_eta) = 1.0_dp
          else
             box%cc(i, j, i_eta) = 0
          end if

          box%cc(i, j, i_pot) = xy(1) + xy(2)
          eta = box%cc(i, j, i_eta)
          pot = box%cc(i, j, i_pot)
          box%cc(i, j, i_rhs) = -4 * eta * pot * inv_dr_sq
       end do
    end do
  end subroutine set_init_cond

  subroutine sides_bc(boxes, id, nb, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, ivs(:)
    integer                     :: nc

    nc = boxes(id)%n_cell

    if (boxes(id)%neighbors(nb) == -1) then
       select case (nb)
       case (a2_nb_lx)
          ! Dirichlet zero
          boxes(id)%cc(0, 1:nc, ivs) = -boxes(id)%cc(1, 1:nc, ivs)
       case (a2_nb_hx)
          ! Dirichlet zero
          boxes(id)%cc(nc+1, 1:nc, ivs) = -boxes(id)%cc(nc, 1:nc, ivs)
       case (a2_nb_ly)
          ! Neumann zero
          boxes(id)%cc(1:nc, 0, ivs) = -boxes(id)%cc(1:nc, 1, ivs)
       case (a2_nb_hy)
          ! Neumann zero
          boxes(id)%cc(1:nc, nc+1, ivs) = -boxes(id)%cc(1:nc, nc, ivs)
       end select
    end if
  end subroutine sides_bc

  subroutine elec_lpl(box, i_in, i_out, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_in, i_out
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq, eta, pot

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2

    do j = 1, nc
       do i = 1, nc
          eta = box%cc(i, j, i_eta)
          pot = box%cc(i, j, i_pot)
          box%cc(i, j, i_out) = inv_dr_sq * ((1-eta) * &
               (box%cc(i-1, j, i_in) + box%cc(i+1, j, i_in) + &
               box%cc(i, j-1, i_in) + box%cc(i, j+1, i_in)) - &
               4 * box%cc(i, j, i_in))
       end do
    end do
  end subroutine elec_lpl

  subroutine elec_gsrb(box, i_phi, i_rhs, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, di(2)
    real(dp)                    :: dxdy, eta

    dxdy = box%dr**2
    nc   = box%n_cell

    ! The parity of redblack_cntr determines which cells we use
    di(1) = iand(redblack_cntr, 1)
    di(2) = ieor(di(1), 1)

    do j = 1, nc, 2
       do i = 1+di(1), nc, 2
          eta = box%cc(i, j, i_eta)
          box%cc(i, j, i_phi) = 0.25_dp * ((1-eta) * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi)) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do

    do j = 2, nc, 2
       do i = 1+di(2), nc, 2
          eta = box%cc(i, j, i_eta)
          box%cc(i, j, i_phi) = 0.25_dp * ((1-eta) * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi)) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine elec_gsrb

end program