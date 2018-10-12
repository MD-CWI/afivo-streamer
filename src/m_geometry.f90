!> Module that provides routines for geometric operations and calculations.
!> Methods and types have a prefix GM_, short for 'geometry'.
module m_geometry

! TODO: Describe methods. Till now those GM methods are not used elsewhere
  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: GM_dist_line
  public :: GM_dist_vec_line
  public :: GM_density_line
  public :: GM_sigmoid
  public :: GM_gaussian
  public :: GM_smoothstep

contains

  !> Compute distance vector between point and its projection onto a line
  !> between r0 and r1
  subroutine GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    integer, intent(in)   :: n_dim
    real(dp), intent(in)  :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp), intent(out) :: dist_vec(n_dim)
    real(dp), intent(out) :: frac !< Fraction [0,1] along line
    real(dp)              :: line_len2

    line_len2 = sum((r1 - r0)**2)
    frac = sum((r - r0) * (r1 - r0))

    if (frac <= 0.0_dp) then
       frac = 0.0_dp
       dist_vec = r - r0
    else if (frac >= line_len2) then
       frac = 1.0_dp
       dist_vec = r - r1
    else
       dist_vec = r - (r0 + frac/line_len2 * (r1 - r0))
       frac = sqrt(frac / line_len2)
    end if
  end subroutine GM_dist_vec_line

  function GM_dist_line(r, r0, r1, n_dim) result(dist)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp)             :: dist, dist_vec(n_dim), frac
    call GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    dist = norm2(dist_vec)
  end function GM_dist_line

  function GM_density_line(r, r0, r1, n_0, n_1, n_dim, width, falloff_t) result(val)
    integer, intent(in)          :: n_dim
    real(dp), intent(in)         :: r(n_dim), r0(n_dim), r1(n_dim), width
    real(dp), intent(in)         :: n_0, n_1
    character(len=*), intent(in) :: falloff_t
    real(dp)                     :: dist, val, dist_vec(n_dim), frac

    call GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    dist = norm2(dist_vec)

    select case (falloff_t)
    case ('sigmoid')
       val  = GM_sigmoid(dist, width)
    case ('gaussian')
       val  = GM_gaussian(dist, width)
    case ('smoothstep')
       val  = GM_smoothstep(dist, width)
    case ('step')
       val  = GM_step(dist, width)
    case ('laser')
       val  = GM_laser(dist_vec, width, n_dim)
    case default
       print *, "GM_density_line: unknown fall-off type: ", trim(falloff_t)
       print *, "Valid options: sigmoid, gaussian, smoothstep, step, laser"
       stop
    end select

    ! Interpolate density between start and endpoint
    val = val * (frac * n_0 + (1-frac) * n_1)
  end function GM_density_line

  function GM_sigmoid(dist, width) result(val)
    real(dp), intent(in) :: dist, width
    real(dp)             :: val, tmp

    tmp = dist / width
    if (tmp > log(0.5_dp * huge(1.0_dp))) then
       val = 0
    else
       val = 2 / (1 + exp(tmp))
    end if
  end function GM_sigmoid

  function GM_gaussian(dist, width) result(val)
    real(dp), intent(in) :: dist, width
    real(dp)             :: val
    val = exp(-(dist/width)**2)
  end function GM_gaussian

  function GM_smoothstep(dist, width) result(val)
    real(dp), intent(in) :: dist, width
    real(dp)             :: val, temp
    if (dist < width) then
       val = 1
    else if (dist < 2 * width) then
       temp = dist/width - 1
       val = (1- (3 * temp**2 - 2 * temp**3))
    else
       val = 0.0_dp
    end if
  end function GM_smoothstep

  function GM_step(dist, width) result(val)
    real(dp), intent(in) :: dist, width
    real(dp)             :: val
    if (dist < width) then
       val = 1
    else
       val = 0.0_dp
    end if
  end function GM_step

  function GM_laser(dist_vec, width, n_dim) result(val)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: dist_vec(n_dim), width
    real(dp)             :: val, xz(2), dy, dxz

    xz(1) = dist_vec(1)
    xz(2) = dist_vec(3)
    dy = abs(dist_vec(2))
    dxz = norm2(xz)

    if (dy < width .and. dxz < width) then
       val = 1
    else
       val = exp(1-(dy**2 + dxz**2)/width**2)
    end if
  end function GM_laser

end module m_geometry
