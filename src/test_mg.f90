program test_mg
  use m_afivo

  implicit none

  type(a2_t)         :: tree
  integer            :: i, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 4
  integer            :: n_boxes_max = 10*1000
  real(dp)           :: dr(2)
  character(len=40)  :: fname, var_names(4)

  type mg2_t
     integer :: i_phi
     integer :: i_phi_old
     integer :: i_rhs
     integer :: i_res
     integer :: n_cycle_up
     integer :: n_cycle_down
  end type mg2_t

  type(mg2_t) :: mg
  mg%i_phi        = 1
  mg%i_phi_old    = 2
  mg%i_res        = 3
  mg%i_rhs        = 4
  mg%n_cycle_up   = 2
  mg%n_cycle_down = 2

  var_names(1) = "phi"
  var_names(2) = "phi_old"
  var_names(3) = "res"
  var_names(4) = "rhs"

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_var_cell=4, n_var_face=0, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  ! Set up the initial conditions
  call a2_set_base(tree, ix_list, nb_list)
  do i = 1, 3
     call a2_adjust_refinement(tree, ref_func_init)
  end do

  call a2_loop_box(tree, set_init_cond)
  call a2_gc_sides(tree, [1,4], mg2d_gc_sides)
  call a2_gc_corners(tree, [1,4], mg2d_gc_corners)

  do i = 1, 10
     write(fname, "(A,I0,A)") "test_mg_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), var_names, i, 0.0_dp)
     call fas_v_cycle(tree, mg, tree%n_levels)
  end do

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    if (box%lvl < 3 .or. (norm2(a2_r_min(box)-2) < 0.5)) then
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
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 1) then
             box%cc(i, j, 4) = 1
          else if (norm2(xy - 2) < 1.1_dp) then
             box%cc(i, j, 4) = (1.1_dp - norm2(xy - 2)) * 10
          else
             box%cc(i, j, 4) = 0
          end if
          box%cc(i, j, 1:3) = 0
       end do
    end do
  end subroutine set_init_cond

  ! On entrance, need valid ghost cell data. On exit, leave valid ghost cell data
  recursive subroutine fas_v_cycle(tree, mg, lvl)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in) :: mg
    integer, intent(in) :: lvl
    integer :: i

    print *, lvl, "pre-relax"
    do i = 1, mg%n_cycle_up
       call gsrb(tree, lvl, mg%i_phi, mg%i_rhs)
    end do

    if (lvl > 1) then
       print *, lvl, "adjust coarse grid"
       call restrict(tree, lvl-1, [mg%i_phi])
       call fill_gc(tree, lvl-1, [mg%i_phi])
       call laplacian(tree, lvl-1, mg%i_rhs, mg%i_phi, .false.)

       call residual(tree, lvl, mg%i_res, mg%i_phi, mg%i_rhs)
       call restrict(tree, lvl-1, [mg%i_res])

       ! rhs_c = laplacian(phi_c) + coarsen(res)
       call add_vars(tree, lvl-1, mg%i_rhs, 1.0_dp, mg%i_res, .false.)

       print *, lvl, "call fas_v_cycle for lvl-1"
       call copy_var(tree, lvl-1, mg%i_phi, mg%i_phi_old)
       call fas_v_cycle(tree, mg, lvl-1)

       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct(tree, lvl, mg%i_phi, mg%i_phi_old)
    end if

    print *, lvl, "post-relax"
    do i = 1, mg%n_cycle_down
       call gsrb(tree, lvl, mg%i_phi, mg%i_rhs)
    end do
  end subroutine fas_v_cycle

  subroutine add_vars(tree, lvl, i_a, a, i_b, use_leafs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: lvl, i_a, i_b
    real(dp), intent(in) :: a
    logical, intent(in) :: use_leafs
    integer                   :: i, id

    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       if (use_leafs .or. a2_has_children(tree%boxes(id))) then
          tree%boxes(id)%cc(:, :, i_a) = &
               tree%boxes(id)%cc(:, :, i_a) &
               + a * tree%boxes(id)%cc(:, :, i_b)
       end if
    end do
  end subroutine add_vars

  ! phi = phi_old + prolong(phi_coarse - phi_old_coarse)
  subroutine correct(tree, lvl, i_phi, i_phi_old)
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: lvl, i_phi, i_phi_old
    integer                   :: i, id

    ! Correct from data on parent level
    do i = 1, size(tree%levels(lvl-1)%ids)
       id = tree%levels(lvl-1)%ids(i)
       if (a2_has_children(tree%boxes(id))) then
          call correct_from_box(tree%boxes, id, i_phi, i_phi_old)
       end if
    end do

    call fill_gc(tree, lvl, [i_phi])
  end subroutine correct

  subroutine correct_from_box(boxes, id, i_phi, i_phi_old)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i_phi, i_phi_old
    real(dp), parameter         :: f1=1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
    integer                     :: nc, i_c, c_id, ix_offset(2)
    integer                     :: i, j, i_c1, i_c2, j_c1, j_c2

    nc = boxes(id)%cfg%n_cell
    do i_c = 1, 4
       c_id = boxes(id)%children(i_c)

       ! Offset of child w.r.t. parent
       ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1)

       ! In these loops, we calculate the closest coarse index (_c1), and the
       ! one-but-closest (_c2). The fine cell lies in between.
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1

             boxes(c_id)%cc(i, j, i_phi) = boxes(c_id)%cc(i, j, i_phi) &
                  + f9 * boxes(id)%cc(i_c1, j_c1, i_phi) &
                  - f9 * boxes(id)%cc(i_c1, j_c1, i_phi_old) &
                  + f3 * boxes(id)%cc(i_c2, j_c1, i_phi) &
                  - f3 * boxes(id)%cc(i_c2, j_c1, i_phi_old) &
                  + f3 * boxes(id)%cc(i_c1, j_c2, i_phi) &
                  - f3 * boxes(id)%cc(i_c1, j_c2, i_phi_old) &
                  + f1 * boxes(id)%cc(i_c2, j_c2, i_phi) &
                  - f1 * boxes(id)%cc(i_c2, j_c2, i_phi_old)
          end do
       end do
    end do
  end subroutine correct_from_box

  subroutine copy_var(tree, lvl, i_from, i_to)
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: lvl, i_from, i_to
    integer                   :: i, id

    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       tree%boxes(id)%cc(:, :, i_to) = tree%boxes(id)%cc(:, :, i_from)
    end do
  end subroutine copy_var

  subroutine restrict(tree, lvl, v_ixs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: lvl, v_ixs(:)
    integer :: i, id

    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       if (a2_has_children(tree%boxes(id))) &
            call a2_restrict_to(tree%boxes, id, v_ixs)
    end do
  end subroutine restrict

  subroutine fill_gc(tree, lvl, v_ixs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: lvl, v_ixs(:)
    integer :: i, id

    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       call a2_gc_box_sides(tree%boxes, id, v_ixs, mg2d_gc_sides)
    end do
    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       call a2_gc_box_corners(tree%boxes, id, v_ixs, mg2d_gc_corners)
    end do
  end subroutine fill_gc

  subroutine gsrb(tree, lvl, i_phi, i_rhs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: lvl, i_phi, i_rhs
    integer :: i, i_red, id

    ! Loop through all fabs in multifab, smoothing red/blue cells only
    do i_red = 1,2
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)
          call gsrb_box(tree%boxes(id), i_phi, i_rhs, i_red==1)
       end do

       ! Communicate updated boundary cells
       call fill_gc(tree, lvl, [i_phi])
    end do
  end subroutine gsrb

  subroutine gsrb_box(box, i_phi, i_rhs, red)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs
    logical, intent(in)         :: red
    integer                     :: i, j, nc, di(2)
    real(dp)                    :: dxdy

    dxdy = product(a2_dr(box))
    nc = box%cfg%n_cell

    if (red) then
       di = (/1,0/)      ! Red cells have i+j = odd
    else
       di = (/0,1/)      ! Black cells have i+j = even
    end if

    do j = 1, nc, 2
       do i = 1+di(1), nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do

    do j = 2, nc, 2
       do i = 1+di(2), nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine gsrb_box

  subroutine laplacian(tree, lvl, i_lpl, i_phi, use_leafs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: lvl, i_lpl, i_phi
    logical, intent(in) :: use_leafs
    integer                   :: i, id

    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       if (use_leafs .or. a2_has_children(tree%boxes(id))) &
            call laplacian_box(tree%boxes(id), i_lpl, i_phi)
    end do
  end subroutine laplacian

  subroutine laplacian_box(box, i_lpl, i_phi)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_lpl, i_phi
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq(2)

    nc = box%cfg%n_cell
    inv_dr_sq = 1 / a2_dr(box)**2

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_lpl) = inv_dr_sq(1) * (box%cc(i-1, j, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i+1, j, i_phi)) + &
               inv_dr_sq(2) * (box%cc(i, j-1, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i, j+1, i_phi))
       end do
    end do
  end subroutine laplacian_box

  subroutine residual(tree, lvl, i_res, i_phi, i_rhs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: lvl, i_res, i_phi, i_rhs
    integer                   :: i, id

    do i = 1, size(tree%levels(lvl)%ids)
       id = tree%levels(lvl)%ids(i)
       call residual_box(tree%boxes(id), i_res, i_phi, i_rhs)
    end do
  end subroutine residual

  subroutine residual_box(box, i_res, i_phi, i_rhs)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_res, i_phi, i_rhs
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq(2)

    nc = box%cfg%n_cell
    inv_dr_sq = 1 / a2_dr(box)**2

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_res) = box%cc(i, j, i_rhs) - ( &
               inv_dr_sq(1) * (box%cc(i-1, j, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i+1, j, i_phi)) + &
               inv_dr_sq(2) * (box%cc(i, j-1, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i, j+1, i_phi)))
          ! print *, box%lvl, i, j, box%cc(i, j, i_res)
       end do
    end do
  end subroutine residual_box

  subroutine mg2d_gc_sides(boxes, id, nb, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, ivs(:)
    integer                     :: n, nc, nv

    nc = boxes(id)%cfg%n_cell
    nv = boxes(id)%cfg%n_var_cell

    select case (nb)
    case (nb_lx)
       if (boxes(id)%neighbors(nb) == -1) then
          boxes(id)%cc(0, 1:nc, ivs) = -boxes(id)%cc(1, 1:nc, ivs)
       else
          call a2_prolong0_to(boxes, id, [0], [(n, n=1,nc)], ivs)
          boxes(id)%cc(0, 1:nc, ivs) = 0.5_dp * boxes(id)%cc(0, 1:nc, ivs) &
               + 0.75_dp * boxes(id)%cc(1, 1:nc, ivs) &
               - 0.25_dp * boxes(id)%cc(2, 1:nc, ivs)
       end if
    case (nb_hx)
       if (boxes(id)%neighbors(nb) == -1) then
          boxes(id)%cc(nc+1, 1:nc, ivs) = -boxes(id)%cc(nc, 1:nc, ivs)
       else
          call a2_prolong0_to(boxes, id, [nc+1], [(n, n=1,nc)], ivs)
          boxes(id)%cc(nc+1, 1:nc, ivs) = 0.5_dp * boxes(id)%cc(nc+1, 1:nc, ivs) &
               + 0.75_dp * boxes(id)%cc(nc, 1:nc, ivs) &
               - 0.25_dp * boxes(id)%cc(nc-1, 1:nc, ivs)
       end if
    case (nb_ly)
       if (boxes(id)%neighbors(nb) == -1) then
          boxes(id)%cc(1:nc, 0, ivs) = -boxes(id)%cc(1:nc, 1, ivs)
       else
          call a2_prolong0_to(boxes, id, [(n, n=1,nc)], [0], ivs)
          boxes(id)%cc(1:nc, 0, ivs) = 0.5_dp * boxes(id)%cc(1:nc, 0, ivs) &
               + 0.75_dp * boxes(id)%cc(1:nc, 1, ivs) &
               - 0.25_dp * boxes(id)%cc(1:nc, 2, ivs)
       end if
    case (nb_hy)
       if (boxes(id)%neighbors(nb) == -1) then
          boxes(id)%cc(1:nc, nc+1, ivs) = -boxes(id)%cc(1:nc, nc, ivs)
       else
          call a2_prolong0_to(boxes, id, [(n, n=1,nc)], [nc+1], ivs)
          boxes(id)%cc(1:nc, nc+1, ivs) = 0.5_dp * boxes(id)%cc(1:nc, nc+1, ivs) &
               + 0.75_dp * boxes(id)%cc(1:nc, nc, ivs) &
               - 0.25_dp * boxes(id)%cc(1:nc, nc-1, ivs)
       end if
    end select
  end subroutine mg2d_gc_sides

  subroutine mg2d_gc_corners(boxes, id, cn, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, cn, ivs(:)
    integer                     :: nc, nv

    nc = boxes(id)%cfg%n_cell
    nv = boxes(id)%cfg%n_var_cell

    select case (cn)
    case (a2_cn_lxly)
       boxes(id)%cc(0, 0, ivs) = boxes(id)%cc(1, 0, ivs) &
            - 0.5_dp * boxes(id)%cc(2, 0, ivs) &
            + boxes(id)%cc(0, 1, ivs) &
            - 0.5_dp * boxes(id)%cc(0, 2, ivs)
    case (a2_cn_hxly)
       boxes(id)%cc(nc+1, 0, ivs) = boxes(id)%cc(nc, 0, ivs) &
            - 0.5_dp * boxes(id)%cc(nc-1, 0, ivs) &
            + boxes(id)%cc(nc+1, 1, ivs) &
            - 0.5_dp * boxes(id)%cc(nc+1, 2, ivs)
    case (a2_cn_lxhy)
       boxes(id)%cc(0, nc+1, ivs) = boxes(id)%cc(0, nc, ivs) &
            - 0.5_dp * boxes(id)%cc(0, nc-1, ivs) &
            + boxes(id)%cc(1, nc+1, ivs) &
            - 0.5_dp * boxes(id)%cc(2, nc+1, ivs)
    case (a2_cn_hxhy)
       boxes(id)%cc(nc+1, nc+1, ivs) = boxes(id)%cc(nc, nc+1, ivs) &
            - 0.5_dp * boxes(id)%cc(nc-1, nc+1, ivs) &
            + boxes(id)%cc(nc+1, nc, ivs) &
            - 0.5_dp * boxes(id)%cc(nc+1, nc-1, ivs)
    end select
  end subroutine mg2d_gc_corners

end program test_mg