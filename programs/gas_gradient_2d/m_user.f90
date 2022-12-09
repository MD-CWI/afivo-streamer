#include "../afivo/src/cpp_macros.h"
!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods
  use m_gas
  use m_streamer

  implicit none
  private

  ! Public methods
  public :: user_initialize

  real(dp)          :: density_ratio = 0.5_dp
  real(dp)          :: shock_width = 0.02_dp
  real(dp)          :: line_coeff(NDIM+1) = 0.0_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_gas_density => gas_density

    call CFG_add_get(cfg, "density_ratio", density_ratio, &
         "Density ratio (<= 1)")
    call CFG_add_get(cfg, "shock_width", shock_width, &
         "Shock width (relative to domain size)")
    call CFG_add_get(cfg, "line_coeff", line_coeff, &
         "Coefficients a, b, c of a line a + bx + cy = 0")

  end subroutine user_initialize

  real(dp) function gas_density(box, IJK)
    use m_gas
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    real(dp)                :: r_rel(NDIM), q, tmp

    ! Get relative coordinates in the range [0, 1]
    r_rel = (af_r_cc(box, [IJK]) - ST_domain_origin)/ ST_domain_len
    q = line_coeff(1) + sum(line_coeff(2:) * r_rel)

    if (q < -shock_width) then
       gas_density = gas_number_density
    else if (q > shock_width) then
       gas_density = gas_number_density * density_ratio
    else
       ! Linear interpolation
       tmp = (q + shock_width) / (2 * shock_width)
       gas_density = gas_number_density * (1 + (density_ratio-1) * tmp)
    end if
  end function gas_density

end module m_user
