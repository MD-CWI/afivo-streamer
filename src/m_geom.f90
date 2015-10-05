module m_geom

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: GM_sigmoid_t    = 1
  integer, parameter :: GM_gaussian_t   = 2
  integer, parameter :: GM_smoothstep_t = 3
  integer, parameter :: GM_step_t       = 4

  ! Public types
  public :: GM_sigmoid_t
  public :: GM_gaussian_t
  public :: GM_smoothstep_t
  public :: GM_step_t

  ! Public methods
  public :: GM_dist_line
  public :: GM_dens_line
  public :: GM_sigmoid
  public :: GM_gaussian
  public :: GM_smoothstep

contains

  !> Compute distance to line between r0 and r1
  function GM_dist_line(r, r0, r1, ndim) result(dist)
    integer, intent(in)  :: ndim
    real(dp), intent(in) :: r(ndim), r0(ndim), r1(ndim)
    real(dp)             :: line_len2, temp, dist
    real(dp)             :: projection(ndim)

    line_len2 = sum((r1 - r0)**2)
    temp = sum((r - r0) * (r1 - r0))

    if (temp <= 0.0_dp) then
       dist = norm2(r-r0)
    else if (temp >= line_len2) then
       dist = norm2(r-r1)
    else
       projection = r0 + temp/line_len2 * (r1 - r0)
       dist = norm2(r-projection)
    end if
  end function GM_dist_line

  function GM_dens_line(r, r0, r1, ndim, width, falloff_t) result(val)
    integer, intent(in)  :: ndim
    real(dp), intent(in) :: r(ndim), r0(ndim), r1(ndim), width
    integer, intent(in)  :: falloff_t
    real(dp)             :: dist, val

    dist = GM_dist_line(r, r0, r1, ndim)

    select case (falloff_t)
    case (GM_sigmoid_t)
       val  = GM_sigmoid(dist, width)
    case (GM_gaussian_t)
       val  = GM_gaussian(dist, width)
    case (GM_smoothstep_t)
       val  = GM_smoothstep(dist, width)
    case (GM_step_t)
       val  = GM_step(dist, width)
    case default
       val  = 0
    end select
  end function GM_dens_line

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

end module m_geom
