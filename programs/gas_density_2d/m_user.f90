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

  character(len=32) :: density_profile = "homogeneous"
  real(dp)          :: density_factor  = 2.0_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_gas_density => gas_density

    call CFG_add_get(cfg, "density_profile", density_profile, &
         "Name of the gas number density profile")
    call CFG_add_get(cfg, "density_factor", density_factor, &
         "Factor for increasing the gas number density")

  end subroutine user_initialize

  real(dp) function gas_density(box, IJK)
    use m_gas
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    real(dp)                :: rz_rel(NDIM), r_rel, z_rel

    ! Get relative r, z coordinates in the range [0, 1]
    rz_rel = (af_r_cc(box, [IJK]) - ST_domain_origin)/ ST_domain_len
    r_rel = rz_rel(1)
    z_rel = rz_rel(2)

    select case (density_profile)
    case ("homogeneous")
       gas_density = gas_number_density
    case ("linear_increase_z")
       ! Linearly increase density with z
       gas_density = gas_number_density * (1 + (density_factor-1) * z_rel)
    case ("linear_decrease_z")
       ! Start from high density and linearly decrease with z
       gas_density = gas_number_density * (1 + (density_factor-1) * (1 - z_rel))
    case default
       error stop "Unknown density_profile specified"
    end select
  end function gas_density

end module m_user
