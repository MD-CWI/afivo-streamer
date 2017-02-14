!> This module can be used to construct solutions consisting of one or more
!> Gaussians.
module m_gaussians

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> A type to store a collection of gaussians in
  type gauss_t
     integer               :: n_gauss  !< Number of gaussians
     integer               :: n_dim    !< Dimensionality
     real(dp), allocatable :: ampl(:)  !< Amplitudes
     real(dp), allocatable :: sigma(:) !< Widths
     real(dp), allocatable :: r0(:,:)  !< Centers
  end type gauss_t

  public :: gauss_t
  public :: gauss_init
  public :: gauss_value
  public :: gauss_gradient
  public :: gauss_laplacian
  public :: gauss_laplacian_cyl
  public :: gauss_4th

contains

  !> Initialize a structure with parameters
  subroutine gauss_init(gs, amplitudes, sigmas, locations)
    type(gauss_t), intent(inout) :: gs             !< Type storing the gaussians
    real(dp), intent(in)         :: amplitudes(:)  !< Their amplitudes
    real(dp), intent(in)         :: sigmas(:)      !< Their widths
    real(dp), intent(in)         :: locations(:,:) !< Their locations

    if (size(locations, 2) /= size(amplitudes) .or. &
         size(sigmas) /= size(amplitudes)) then
       stop "gauss_init: arguments do not match in size"
    end if

    gs%n_gauss = size(amplitudes)
    gs%n_dim = size(locations, 1)

    allocate(gs%ampl(gs%n_gauss))
    allocate(gs%sigma(gs%n_gauss))
    allocate(gs%r0(gs%n_dim, gs%n_gauss))

    gs%ampl = amplitudes
    gs%sigma = sigmas
    gs%r0 = locations
  end subroutine gauss_init

  !> Return the value of the sum of gaussians at r
  real(dp) function gauss_value(gs, r)
    type(gauss_t), intent(in) :: gs
    real(dp), intent(in)      :: r(gs%n_dim)
    integer                   :: n

    gauss_value = 0
    do n = 1, gs%n_gauss
       gauss_value = gauss_value + gauss_single(gs, r, n)
    end do
  end function gauss_value

  !> Return the value of a single gaussian at r
  real(dp) function gauss_single(gs, r, ix)
    type(gauss_t), intent(in) :: gs
    real(dp), intent(in)      :: r(gs%n_dim)
    integer, intent(in)       :: ix
    real(dp)                  :: xrel(gs%n_dim)

    xrel = (r-gs%r0(:, ix)) / gs%sigma(ix)
    gauss_single = exp(-sum(xrel**2))
  end function gauss_single

  subroutine gauss_gradient(gs, r, grad)
    type(gauss_t), intent(in) :: gs
    real(dp), intent(in)      :: r(gs%n_dim)
    real(dp), intent(out)     :: grad(gs%n_dim)
    integer                   :: ix
    real(dp)                  :: xrel(gs%n_dim)

    grad = 0
    do ix = 1, gs%n_gauss
       xrel = (r-gs%r0(:, ix)) / gs%sigma(ix)
       grad = grad - 2 * xrel/gs%sigma(ix) * &
            gauss_single(gs, r, ix)
    end do
  end subroutine gauss_gradient

  !> Summed Laplacian of the gaussians in Cartesian coordinates
  real(dp) function gauss_laplacian(gs, r)
    type(gauss_t), intent(in) :: gs
    real(dp), intent(in)      :: r(gs%n_dim)
    integer                   :: ix
    real(dp)                  :: xrel(gs%n_dim)

    gauss_laplacian = 0
    do ix = 1, gs%n_gauss
       xrel = (r-gs%r0(:, ix)) / gs%sigma(ix)
       gauss_laplacian = gauss_laplacian + 4/gs%sigma(ix)**2 * &
            (sum(xrel**2) - 0.5_dp * gs%n_dim) * gauss_single(gs, r, ix)
    end do
  end function gauss_laplacian

  !> Summed Laplacian of the gaussians in (r,z) coordinates
  real(dp) function gauss_laplacian_cyl(gs, r)
    type(gauss_t), intent(in) :: gs
    real(dp), intent(in)      :: r(gs%n_dim)
    integer :: ix
    real(dp)                  :: xrel(gs%n_dim)

    gauss_laplacian_cyl = 0
    do ix = 1, gs%n_gauss
       xrel = (r-gs%r0(:, ix)) / gs%sigma(ix)
       gauss_laplacian_cyl = gauss_laplacian_cyl + 4/gs%sigma(ix)**2 * &
            (sum(xrel**2) - 1 - 0.5_dp * (r(1)-gs%r0(1, ix))/r(1)) * &
            gauss_single(gs, r, ix)
    end do
  end function gauss_laplacian_cyl

  !> Fourth derivative of the gaussians in Cartesian coordinates
  real(dp) function gauss_4th(gs, r)
    type(gauss_t), intent(in) :: gs
    real(dp), intent(in)      :: r(gs%n_dim)
    integer                   :: ix
    real(dp)                  :: xrel(gs%n_dim), d4(gs%n_dim)

    d4 = 0
    do ix = 1, gs%n_gauss
       xrel = (r-gs%r0(:, ix)) / gs%sigma(ix)
       d4 = d4 + gauss_single(gs, r, ix) / gs%sigma(ix)**4 * &
            (16 * xrel**4 - 48 * xrel**2  + 12)
    end do
    gauss_4th = norm2(d4)
  end function gauss_4th

end module m_gaussians
