module m_mg2d
  use m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type mg2_t
     private
     integer  :: i_phi
     integer  :: i_phi_old
     integer  :: i_rhs
     integer  :: i_res
     integer  :: n_cycle_down
     integer  :: n_cycle_up
     real(dp) :: base_res_reduction
     procedure(a2_subr_gc), pointer, nopass :: sides_bc, corners_bc
  end type mg2_t

  ! Public types
  public :: mg2_t

  ! Public methods
  public :: mg2d_set
  public :: mg2d_restrict_rhs
  public :: mg2d_fas_fmg
  public :: mg2d_fas_vcycle

contains

  subroutine mg2d_set(mg, i_phi, i_phi_old, i_rhs, i_res, &
       n_cycle_down, n_cycle_up, sides_bc, corners_bc)
    type(mg2_t), intent(out) :: mg
    integer, intent(in)      :: i_phi, i_phi_old, i_rhs, i_res
    integer, intent(in)      :: n_cycle_down, n_cycle_up
    procedure(a2_subr_gc)    :: sides_bc, corners_bc
    mg%i_phi        = i_phi
    mg%i_phi_old    = i_phi_old
    mg%i_rhs        = i_rhs
    mg%i_res        = i_res
    mg%n_cycle_down = n_cycle_down
    mg%n_cycle_up   = n_cycle_up
    mg%sides_bc     => sides_bc
    mg%corners_bc   => corners_bc
  end subroutine mg2d_set

  subroutine mg2d_restrict_rhs(tree, mg)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(inout) :: mg
    integer :: lvl, i, id

    do lvl = tree%n_lvls-1, 1, -1
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_restrict_to(tree%boxes, id, [mg%i_rhs])
       end do
    end do
  end subroutine mg2d_restrict_rhs

  ! Need valid ghost cells on input, has valid gc on output
  subroutine mg2d_fas_fmg(tree, mg)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in)   :: mg
    integer                   :: i, id, lvl

    do lvl = 1, tree%n_lvls
       ! Store phi in phi_old
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          tree%boxes(id)%cc(:, :, mg%i_phi_old) = &
               tree%boxes(id)%cc(:, :, mg%i_phi)
       end do

       if (lvl > 1) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          do i = 1, size(tree%lvls(lvl-1)%parents)
             id = tree%lvls(lvl-1)%parents(i)
             call correct_from_box(tree%boxes, id, mg%i_phi, mg%i_phi_old)
          end do

          ! Update ghost cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       endif

       call mg2d_fas_vcycle(tree, mg, lvl) ! Perform V-cycle
    end do
  end subroutine mg2d_fas_fmg


  ! Modified from gmres code found at:
  ! http://people.sc.fsu.edu/~jburkardt/f_src/mgmres/mgmres.html
  ! Which has the GNU LGPL license.
  ! Original C version by Lili Ju., FORTRAN90 version by John Burkardt.
  ! Modifications: Jannis Teunissen
  subroutine mg2d_gmres(x, rhs, a_times_x, tree, mg, max_outer, max_inner, &
       tol_abs, tol_rel)
    integer, intent(in)       :: max_outer, max_inner
    real(dp), intent(in)      :: tol_abs, tol_rel, rhs(:)
    real(dp), intent(inout)   :: x(:)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in)   :: mg

    interface
       subroutine a_times_x(x, ax, tree, mg)
         import
         real(dp), intent(in)      :: x(:)
         real(dp), intent(out)     :: ax(:)
         type(a2_t), intent(inout) :: tree
         type(mg2_t), intent(in)   :: mg
       end subroutine a_times_x
    end interface

    real(dp), parameter :: delta   = 1.0e-3_dp
    logical, parameter  :: verbose = .true.
    logical             :: finished
    integer             :: i, j, k, itr, itr_used
    real(dp)            :: prev_norm, r_tmp, rho, rho_tol, radius
    real(dp)            :: res(size(x))
    real(dp)            :: gcos(max_inner)
    real(dp)            :: gsin(max_inner)
    real(dp)            :: g(max_inner+1)
    real(dp)            :: y(max_inner+1)
    real(dp)            :: h(max_inner+1, max_inner)
    real(dp)            :: v(size(x), max_inner+1)

    itr_used = 0

    do itr = 1, max_outer
       call a_times_x(x, res, tree, mg)

       res      = rhs - res
       rho      = norm2(res)

       ! Use first residual to set tolerance
       if (itr == 1) rho_tol = rho * tol_rel

       ! Check whether we can stop already
       if (rho <= rho_tol .and. rho <= tol_abs) exit

       r_tmp  = 1/rho
       v(:,1) = res * r_tmp
       g(1)   = rho
       g(2:)  = 0.0_dp
       h(:,:) = 0.0_dp

       if (verbose) print *, "outer", itr, "norm residual:", rho

       finished = .false.
       k        = 0

       do
          if (finished) exit
          k = k + 1

          call a_times_x(v(:,k), v(:,k+1), tree, mg)
          prev_norm = norm2(v(:,k+1))

          ! Orthogonalize new vector
          do j = 1, k
             h(j,k)   = dot_product(v(:,k+1), v(:,j))
             v(:,k+1) = v(:,k+1) - h(j,k) * v(:,j)
          end do

          ! Store norm of new vector
          h(k+1,k) = norm2(v(:,k+1))

          ! If the orthogonalized vector norm is very small compared to the
          ! initial norm, orthogonalize again to improve accuracy
          if (prev_norm + delta * h(k+1,k) == prev_norm) then
             do j = 1, k
                r_tmp    = dot_product(v(:,k+1), v(:,j))
                h(j,k)   = h(j,k) + r_tmp
                v(:,k+1) = v(:,k+1) - r_tmp * v(:,j)
             end do
             h(k+1,k) = norm2(v(:,k+1))
          end if

          ! Normalize new vector, but avoid division by zero. If division by
          ! zero would occur, we will exit at the next convergence check.
          if (h(k+1, k) > epsilon(1.0_dp)) then
             r_tmp = 1/h(k+1,k)
             v(:,k+1) = v(:,k+1) * r_tmp
          else
             finished = .true.
          end if

          if (k > 1) then
             y(1:k+1) = h(1:k+1,k)

             do j = 1, k-1
                call mult_givens(gcos(j), gsin(j), j, y(1:k+1))
             end do

             h(1:k+1,k) = y(1:k+1)
          end if

          ! Compute givens rotation angle cos/sin
          radius   = hypot(h(k,k), h(k+1,k))
          r_tmp    = 1/radius
          gcos(k)  = h(k,k) * r_tmp
          gsin(k)  = -h(k+1,k) * r_tmp
          h(k,k)   = radius
          h(k+1,k) = 0.0_dp

          call mult_givens(gcos(k), gsin(k), k, g(1:k+1))

          rho      = abs(g(k+1))
          finished = finished .or. (k == max_inner) .or. &
               (rho <= rho_tol .and. rho <= tol_abs)

          if (verbose) print *, "inner", k, "norm residual:", rho
       end do

       ! Update solution guess x
       y(k) = g(k) / h(k,k)

       do i = k-1, 1, -1
          y(i) = (g(i) - dot_product(h(i,i+1:k), y(i+1:k))) / h(i,i)
       end do

       do i = 1, size(x)
          x(i) = x(i) + dot_product(v(i,1:k), y(1:k))
       end do
    end do
  end subroutine mg2d_gmres

  ! Perform givens rotation
  subroutine mult_givens(gcos, gsin, k, my_vec)
    integer, intent(in)     :: k
    real(dp), intent(in)    :: gcos, gsin
    real(dp), intent(inout) :: my_vec(:)
    real(dp)                :: g1, g2

    g1          = gcos * my_vec(k) - gsin * my_vec(k+1)
    g2          = gsin * my_vec(k) + gcos * my_vec(k+1)
    my_vec(k)   = g1
    my_vec(k+1) = g2
  end subroutine mult_givens

  subroutine mg2d_lpl_vec(phi, lpl, tree, mg)
    real(dp), intent(in) :: phi(:)
    real(dp), intent(out) :: lpl(:)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in)   :: mg
    integer                   :: i, id, nc, nc2

    nc     = tree%cfg%n_cell
    nc2    = nc**2

    ! Put the vector phi back in the boxes
    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       tree%boxes(id)%cc(1:nc, 1:nc, mg%i_phi) = &
            reshape(phi((i-1)*nc2+1:i*nc2), [nc,nc])
    end do

    ! Communicate updated boundary cells
    call mg2d_fill_gc(tree%boxes, tree%lvls(1)%ids, [mg%i_phi], &
         mg%sides_bc, mg%corners_bc)

    ! Calculate laplacian
    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       call laplacian_box(tree%boxes(id), mg%i_rhs, mg%i_phi)
    end do

    ! Put the result in lpl
    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       lpl((i-1)*nc2+1:i*nc2) = &
            reshape(tree%boxes(id)%cc(1:nc, 1:nc, mg%i_phi), [nc2])
    end do

  end subroutine mg2d_lpl_vec

  ! On entrance, need valid ghost cell data. On exit, leave valid ghost cell data
  recursive subroutine mg2d_fas_vcycle(tree, mg, lvl)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in)   :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, n, n_base, nc, nc2
    real(dp), allocatable :: phi_base(:), rhs_base(:)

    if (lvl == 1) then
       ! Use gmres to solve base level
       print *, "Solving base level"
       nc     = tree%cfg%n_cell
       nc2    = nc**2
       n_base = size(tree%lvls(lvl)%ids)
       allocate(phi_base(n_base * nc2))
       allocate(rhs_base(n_base * nc2))

       ! Pack phi at the base level into a vector
       do i = 1, n_base
          id = tree%lvls(lvl)%ids(i)
          phi_base((i-1)*nc2+1:i*nc2) = &
               reshape(tree%boxes(id)%cc(1:nc, 1:nc, mg%i_phi), [nc2])
       end do

       call mg2d_gmres(phi_base, rhs_base, mg2d_lpl_vec, tree, mg, 1, nc2, &
            huge(1.0_dp), mg%base_res_reduction)

       deallocate(phi_base)
       deallocate(rhs_base)
    else
       ! Downwards relaxation. Half the cycle are "red", half are "black".
       do n = 1, 2 * mg%n_cycle_down
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call gsrb_box(tree%boxes(id), mg%i_phi, mg%i_rhs, n)
          end do

          ! Communicate updated boundary cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       end do

       ! Calculate residual at current lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call residual_box(tree%boxes(id), mg%i_res, mg%i_phi, mg%i_rhs)
       end do

       ! Restrict phi
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call a2_restrict_to(tree%boxes, id, [mg%i_phi])
       end do

       ! Have to update ghost cells for phi_c (todo: not everywhere?)
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl-1)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call laplacian_box(tree%boxes(id), mg%i_rhs, mg%i_phi)
          call a2_restrict_to_iadd(tree%boxes, id, mg%i_res, mg%i_rhs)
       end do

       ! Store current coarse phi in phi_old
       do i = 1, size(tree%lvls(lvl-1)%ids)
          id = tree%lvls(lvl-1)%ids(i)
          tree%boxes(id)%cc(:, :, mg%i_phi_old) = &
               tree%boxes(id)%cc(:, :, mg%i_phi)
       end do

       call mg2d_fas_vcycle(tree, mg, lvl-1)

       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call correct_from_box(tree%boxes, id, mg%i_phi, mg%i_phi_old)
       end do

       ! Have to fill ghost cells again (todo: not everywhere?)
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Upwards relaxation
       do n = 1, 2 * mg%n_cycle_down
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call gsrb_box(tree%boxes(id), mg%i_phi, mg%i_rhs, n)
          end do

          ! Communicate updated boundary cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       end do
    end if

    ! Calculate residual at current lvl for output
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       call residual_box(tree%boxes(id), mg%i_res, mg%i_phi, mg%i_rhs)
    end do
  end subroutine mg2d_fas_vcycle

  subroutine mg2d_fill_gc(boxes, ids, ivs, sides_bc, corners_bc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), ivs(:)
    procedure(a2_subr_gc)       :: sides_bc, corners_bc
    integer                     :: i

    do i = 1, size(ids)
       call a2_gc_box_sides(boxes, ids(i), ivs, &
            a2_sides_extrap, sides_bc)
    end do
    do i = 1, size(ids)
       call a2_gc_box_corners(boxes, ids(i), ivs, &
            a2_corners_extrap, corners_bc)
    end do
  end subroutine mg2d_fill_gc

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
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
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

             boxes(c_id)%cc(i, j, i_phi) = boxes(c_id)%cc(i, j, i_phi) &
                  + f9 * (boxes(id)%cc(i_c1, j_c1, i_phi) &
                  - boxes(id)%cc(i_c1, j_c1, i_phi_old)) &
                  + f3 * (boxes(id)%cc(i_c2, j_c1, i_phi) &
                  - boxes(id)%cc(i_c2, j_c1, i_phi_old) &
                  + boxes(id)%cc(i_c1, j_c2, i_phi) &
                  - boxes(id)%cc(i_c1, j_c2, i_phi_old)) &
                  + f1 * (boxes(id)%cc(i_c2, j_c2, i_phi) &
                  - boxes(id)%cc(i_c2, j_c2, i_phi_old))
          end do
       end do
    end do
  end subroutine correct_from_box

  subroutine gsrb_box(box, i_phi, i_rhs, redblack_cntr)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    integer                     :: i, j, nc, di(2)
    real(dp)                    :: dxdy

    dxdy = product(a2_dr(box))
    nc   = box%cfg%n_cell

    ! The parity of redblack_cntr determines which cells we use
    di(1) = iand(redblack_cntr, 1)
    di(2) = ieor(di(1), 1)

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
       end do
    end do
  end subroutine residual_box

end module m_mg2d