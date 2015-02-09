!> \example test_mg2_2d.f90
!> Example showing how to use m_mg_2d for a jumping coefficient (e.g.,
!> a jump in dielectric constant in electrostatics)
program test_mg3_2d
  use m_afivo_2d
  use m_mg_2d

  implicit none

  integer, parameter :: dp           = kind(0.0d0)
  integer, parameter :: box_size     = 2
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_res = 4
  integer, parameter :: i_eps = 5

  type(a2_t)         :: tree
  integer            :: i, id, n_changes
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  integer            :: n_boxes_max  = 20*1000
  integer            :: n_lvls_max   = 20
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(5)
  type(mg2_t)        :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_res) = "res"
  var_names(i_eps) = "eps"

  dr = 1.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_lvls_max, n_boxes_max, box_size, n_var_cell=5, &
       n_var_face=0, dr = dr, r_min = [0.0_dp, 0.0_dp], coarsen_to=-1)

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

  ! Set the multigrid options
  mg%i_phi        = i_phi
  mg%i_tmp        = i_tmp
  mg%i_rhs        = i_rhs
  mg%i_res        = i_res
  mg%n_cycle_down = 2
  mg%n_cycle_up   = 2
  mg%n_cycle_base = 2
  mg%sides_bc     => sides_bc
  mg%box_op       => my_lpl_box
  mg%box_gsrb     => my_gsrb_lpl_box
  mg%box_corr     => my_corr_lpl_box

  call mg2_init_mg(mg)

  ! Restrict from children recursively
  call a2_restrict_tree(tree, i_rhs)
  call a2_restrict_tree(tree, i_phi)

  !$omp parallel
  do i = 1, 10
     ! call mg2_fas_vcycle(tree, mg, tree%n_lvls)
     call mg2_fas_fmg(tree, mg)
     write(fname, "(A,I0,A)") "test_mg3_2d_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), var_names, i, 0.0_dp)
  end do
  !$omp end parallel

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    integer                  :: n

    ref_func_init = a5_rm_ref

    if (box%lvl < 6) then
       ref_func_init = a5_do_ref
    else
       ref_func_init = a5_kp_ref
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
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_rhs) = 0.0_dp
          if (all(xy > 0.5) .or. all(xy < 0.5)) then
             box%cc(i, j, i_eps) = 1.0_dp
          else
             box%cc(i, j, i_eps) = 1000.0_dp
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    real(dp)                    :: xy(2)
    integer                     :: n, nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       ! Neumann zero
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       ! Neumann zero
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       ! Dirichlet zero
       boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       ! Dirichlet one
       boxes(id)%cc(1:nc, nc+1, iv) = 2 - boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc

  subroutine my_gsrb_lpl_box(box, i_phi, i_rhs, redblack_cntr)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_phi !< Index of solution variable
    integer, intent(in)         :: i_rhs !< Index of right-hand side
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    integer                     :: i, i0, j, nc
    real(dp)                    :: dx2, u(5), a(5), c(4)

    dx2 = box%dr**2
    nc  = box%n_cell

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          u(1) = box%cc(i, j, i_phi)
          u(2:3) = box%cc(i-1:i+2:2, j, i_phi)
          u(4:5) = box%cc(i, j-1:j+2:2, i_phi)
          a(1) = box%cc(i, j, i_eps)
          a(2:3) = box%cc(i-1:i+2:2, j, i_eps)
          a(4:5) = box%cc(i, j-1:j+2:2, i_eps)
          c(:) = 2 * a(1) * a(2:) / (a(1) + a(2:))

          box%cc(i, j, i_phi) = &
               (sum(c(:) * u(2:)) - dx2 * box%cc(i, j, i_rhs)) / sum(c)
       end do
    end do
  end subroutine my_gsrb_lpl_box

  !> Perform Laplacian operator on a box
  subroutine my_lpl_box(box, i_in, i_out)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_in !< Index of variable to take Laplacian of
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq, u(5), a(5)

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2

    do j = 1, nc
       do i = 1, nc
          u(1) = box%cc(i, j, i_in)
          u(2:3) = box%cc(i-1:i+2:2, j, i_in)
          u(4:5) = box%cc(i, j-1:j+2:2, i_in)
          a(1) = box%cc(i, j, i_eps)
          a(2:3) = box%cc(i-1:i+2:2, j, i_eps)
          a(4:5) = box%cc(i, j-1:j+2:2, i_eps)

          box%cc(i, j, i_out) = inv_dr_sq * 2 * &
               sum(a(1)*a(2:)/(a(1) + a(2:)) * (u(2:) - u(1)))
       end do
    end do
  end subroutine my_lpl_box

  subroutine my_corr_lpl_box(box_p, box_c, i_phi, i_corr)
    type(box2_t), intent(inout) :: box_c
    type(box2_t), intent(in)    :: box_p
    integer, intent(in)         :: i_phi, i_corr
    integer                     :: ix_offset(2)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp) :: u(3), a(3)

    nc = box_c%n_cell
    ix_offset = a2_get_child_offset(box_c)

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          u(1) = box_p%cc(i_c1, j_c1, i_corr)
          u(2) = box_p%cc(i_c2, j_c1, i_corr)
          u(3) = box_p%cc(i_c1, j_c2, i_corr)
          a(1) = box_p%cc(i_c1, j_c1, i_eps)
          a(2) = box_p%cc(i_c2, j_c1, i_eps)
          a(3) = box_p%cc(i_c1, j_c2, i_eps)

          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + 0.5_dp * ( &
               (a(1)*u(1)+a(2)*u(2)) / (a(1) + a(2)) + &
               (a(1)*u(1)+a(3)*u(3)) / (a(1) + a(3)))
       end do
    end do
  end subroutine my_corr_lpl_box

end program test_mg3_2d
