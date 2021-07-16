#include "../afivo/src/cpp_macros.h"
!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods
  use m_streamer

  implicit none
  private

  ! Public methods
  public :: user_initialize

  ! What kind of dielectric to use
  character(len=20) :: dielectric_type = "top"
  real(dp) :: dielectric_eps = 2.0_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_initial_conditions => my_init_cond

    call CFG_add_get(cfg, "dielectric_type", dielectric_type, &
         "What kind of dielectric to use")
    call CFG_add_get(cfg, "dielectric_eps", dielectric_eps, &
         "The dielectric permittivity")

  end subroutine user_initialize

  subroutine my_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    select case (dielectric_type)
       case ("top")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(2) > 0.75_dp * ST_domain_len(2)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("bottom")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(2) < 0.25_dp * ST_domain_len(2)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("top_bottom")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(2) > 0.75_dp * ST_domain_len(2) .or. &
                  rr(2) < 0.25_dp * ST_domain_len(2)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("left")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(1) < 0.25_dp * ST_domain_len(1)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("rod")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(1) < 0.125_dp * ST_domain_len(1)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("hollow_rod")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(1) > 0.0625_dp * ST_domain_len(1) .and. &
                  rr(1) < 0.125_dp * ST_domain_len(1)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("left_right")
          do KJI_DO(0,nc+1)
             rr   = af_r_cc(box, [IJK])
             if (rr(1) < 0.25_dp * ST_domain_len(1) .or. &
                  rr(1) > 0.75_dp * ST_domain_len(1)) then
                box%cc(IJK, i_eps) = dielectric_eps
                box%cc(IJK, i_electron) = 0.0_dp
                box%cc(IJK, i_1pos_ion) = 0.0_dp
             else
                box%cc(IJK, i_eps) = 1.0_dp
             end if
          end do; CLOSE_DO
       case ("gas")
          do KJI_DO(0,nc+1)
                box%cc(IJK, i_eps) = 1.0_dp
          end do; CLOSE_DO
       case default
          error stop "Unknown dielectric_type"
       end select
  end subroutine my_init_cond

end module m_user
