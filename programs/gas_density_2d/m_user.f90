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

  character(len=32) :: density_profile_z = "homogeneous"
  character(len=32) :: density_profile_r = "homogeneous"
  real(dp)          :: z_density_ratio  = 0.0_dp
  real(dp)          :: r_reduction  = 0.5_dp
  real(dp)          :: r_width  = 0.1_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_gas_density => gas_density

    call CFG_add_get(cfg, "density_profile_z", density_profile_z, &
         "Name of the gas number density profile in the z direction")
    call CFG_add_get(cfg, "density_profile_r", density_profile_r, &
         "Name of the gas number density profile in the r direction")
    call CFG_add_get(cfg, "z_density_ratio", z_density_ratio, &
         "Density ratio in the z direction")
    call CFG_add_get(cfg, "r_reduction", r_reduction, &
         "Reduction of the gas number density on the axis")
    call CFG_add_get(cfg, "r_width", r_width, &
         "Width of the profile in the r direction")

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

    select case (density_profile_z)
    case ("homogeneous")
       gas_density = gas_number_density
    case ("linear_z")
       ! Linearly trend in z
       gas_density = gas_number_density * (1 + (z_density_ratio-1) * z_rel) / &
            max(1.0_dp, abs(z_density_ratio))
    case default
       error stop "Unknown density_profile_z specified"
    end select

    select case (density_profile_r)
    case ("homogeneous")
       continue
    case ("gaussian")
       gas_density = gas_density * (1 - r_reduction * exp(-(r_rel/r_width)**2))
    case ("step")
       if (r_rel < r_width) then
          ! Reduce gas density inside channel
          gas_density = r_reduction * gas_density
       end if
    case default
       error stop "Unknown density_profile_r specified"
    end select
  end function gas_density

end module m_user
