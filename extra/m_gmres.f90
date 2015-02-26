module m_gmres

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  abstract interface
     subroutine subr_ax(x, ax)
       import
       real(dp), intent(in) :: x(:)
       real(dp), intent(out) :: ax(:)
     end subroutine subr_ax
  end interface

  ! Public interfac
  public :: subr_ax

  ! Public methods
  public :: gmr_gmres

contains

  ! Modified from gmres code found at:
  ! http://people.sc.fsu.edu/~jburkardt/f_src/mgmres/mgmres.html
  ! Which has the GNU LGPL license.
  ! Original C version by Lili Ju., FORTRAN90 version by John Burkardt.
  ! Modifications: Jannis Teunissen
  subroutine gmr_gmres(x, rhs, a_times_x, max_outer, max_inner, &
       tol_abs, tol_rel)
    integer, intent(in)     :: max_outer, max_inner
    real(dp), intent(in)    :: tol_abs, tol_rel, rhs(:)
    real(dp), intent(inout) :: x(:)
    procedure(subr_ax)      :: a_times_x

    real(dp), parameter     :: delta   = 1.0e-3_dp
    logical, parameter      :: verbose = .false.

    logical                 :: finished
    integer                 :: i, j, k, itr, itr_used

    real(dp)                :: prev_norm, r_tmp, rho, rho_tol, radius
    real(dp)                :: res(size(x))
    real(dp)                :: gcos(max_inner)
    real(dp)                :: gsin(max_inner)
    real(dp)                :: g(max_inner+1)
    real(dp)                :: y(max_inner+1)
    real(dp)                :: h(max_inner+1, max_inner)
    real(dp)                :: v(size(x), max_inner+1)

    itr_used = 0

    do itr = 1, max_outer
       call a_times_x(x, res)

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

          call a_times_x(v(:,k), v(:,k+1))
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
  end subroutine gmr_gmres

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

end module m_gmres