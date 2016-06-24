!> This module can be used to construct solutions consisting of one or more
!> Gaussians.
!>
!> Author: Jannis Teunissen
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
  public :: gauss_laplacian_cylindrical
  public :: gauss_4th

contains

  !> Initialize a structure with parameters
  subroutine gauss_init(gaussian, amplitudes, sigmas, locations)
    type(gauss_t), intent(inout) :: gaussian       !< Type storing the gaussians
    real(dp), intent(in)         :: amplitudes(:)  !< Their amplitudes
    real(dp), intent(in)         :: sigmas(:)      !< Their widths
    real(dp), intent(in)         :: locations(:,:) !< Their locations

    if (size(locations, 2) /= size(amplitudes) .or. &
         size(sigmas) /= size(amplitudes)) then
       stop "gauss_init: arguments do not match in size"
    end if

    gaussian%n_gauss = size(amplitudes)
    gaussian%n_dim = size(locations, 1)

    allocate(gaussian%ampl(gaussian%n_gauss))
    allocate(gaussian%sigma(gaussian%n_gauss))
    allocate(gaussian%r0(gaussian%n_dim, gaussian%n_gauss))

    gaussian%ampl = amplitudes
    gaussian%sigma = sigmas
    gaussian%r0 = locations
  end subroutine gauss_init

  !> Return the value of the sum of gaussians at r
  real(dp) function gauss_value(gaussian, r)
    type(gauss_t), intent(in) :: gaussian
    real(dp), intent(in)      :: r(gaussian%n_dim)
    integer                   :: n

    gauss_value = 0
    do n = 1, gaussian%n_gauss
       gauss_value = gauss_value + gauss_single(gaussian, r, n)
    end do
  end function gauss_value

  !> Return the value of a single gaussian at r
  real(dp) function gauss_single(gaussian, r, ix)
    type(gauss_t), intent(in) :: gaussian
    real(dp), intent(in)      :: r(gaussian%n_dim)
    integer, intent(in)       :: ix
    real(dp)                  :: xrel(gaussian%n_dim)

    xrel = (r-gaussian%r0(:, ix)) / gaussian%sigma(ix)
    gauss_single = exp(-sum(xrel**2))
  end function gauss_single

  subroutine gauss_gradient(gaussian, r, gradient)
    type(gauss_t), intent(in) :: gaussian
    real(dp), intent(in)      :: r(gaussian%n_dim)
    real(dp), intent(out)     :: gradient(gaussian%n_dim)
    integer                   :: ix
    real(dp)                  :: xrel(gaussian%n_dim)

    gradient = 0
    do ix = 1, gaussian%n_gauss
       xrel = (r-gaussian%r0(:, ix)) / gaussian%sigma(ix)
       gradient = gradient - 2 * xrel/gaussian%sigma(ix) * &
            gauss_single(gaussian, r, ix)
    end do
  end subroutine gauss_gradient

  !> Summed Laplacian of the gaussians in Cartesian coordinates
  real(dp) function gauss_laplacian(gaussian, r)
    type(gauss_t), intent(in) :: gaussian
    real(dp), intent(in)      :: r(gaussian%n_dim)
    integer                   :: ix
    real(dp)                  :: xrel(gaussian%n_dim)

    gauss_laplacian = 0
    do ix = 1, gaussian%n_gauss
       xrel = (r-gaussian%r0(:, ix)) / gaussian%sigma(ix)
       gauss_laplacian = gauss_laplacian + 4/gaussian%sigma(ix)**2 * &
            (sum(xrel**2) - 0.5_dp * gaussian%n_dim) * gauss_single(gaussian, r, ix)
    end do
  end function gauss_laplacian

  !> Summed Laplacian of the gaussians in (r,z) coordinates
  real(dp) function gauss_laplacian_cylindrical(gaussian, r)
    type(gauss_t), intent(in) :: gaussian
    real(dp), intent(in)      :: r(gaussian%n_dim)
    integer :: ix
    real(dp)                  :: xrel(gaussian%n_dim)

    gauss_laplacian_cylindrical = 0
    do ix = 1, gaussian%n_gauss
       xrel = (r-gaussian%r0(:, ix)) / gaussian%sigma(ix)
       gauss_laplacian_cylindrical = gauss_laplacian_cylindrical + 4/gaussian%sigma(ix)**2 * &
            (sum(xrel**2) - 1 - 0.5_dp * (r(1)-gaussian%r0(1, ix))/r(1)) * &
            gauss_single(gaussian, r, ix)
    end do
  end function gauss_laplacian_cylindrical

  !> Fourth derivative of the gaussians in Cartesian coordinates
  real(dp) function gauss_4th(gaussian, r)
    type(gauss_t), intent(in) :: gaussian
    real(dp), intent(in)      :: r(gaussian%n_dim)
    integer :: ix
    real(dp)                  :: xrel(gaussian%n_dim)

    gauss_4th = 0
    do ix = 1, gaussian%n_gauss
       xrel = (r-gaussian%r0(:, ix)) / gaussian%sigma(ix)
       gauss_4th = gauss_4th + gauss_single(gaussian, r, ix) / gaussian%sigma(ix)**4 * &
            (16 * sum(xrel**4) - 48 * sum(xrel**2) + gaussian%n_dim * 12)
    end do
  end function gauss_4th

end module m_gaussians
