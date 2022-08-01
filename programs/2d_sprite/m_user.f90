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

  ! Wait-Spies model
  real(dp) :: e_decay_height = 2.86e3_dp

  real(dp) :: scale_height   = 7.2e3_dp
  real(dp) :: n_e0           = 1e4_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_initial_conditions => my_init_cond
    user_gas_density => gas_density

  end subroutine user_initialize

  pure real(dp) function gas_density(box, IJK)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    real(dp)                :: rr(NDIM)

    rr = af_r_cc(box, [IJK])
    gas_density = 2.5e25 * exp(-rr(NDIM) / scale_height)
  end function gas_density

  subroutine my_init_cond(box)
    use m_init_cond
    use m_geometry
    type(box_t), intent(inout) :: box
    integer                    :: i, j, n
    real(dp)                   :: rr(2), n_e, density

    do j = 0, box%n_cell+1
       do i = 0, box%n_cell+1
          rr = af_r_cc(box, [i, j])
          n_e = n_e0 * exp((rr(2) - 60e3) / e_decay_height)

          box%cc(i, j, i_electron) = n_e
          box%cc(i, j, i_1pos_ion) = n_e

          do n = 1, init_conds%n_cond
             density = GM_density_line(rr, init_conds%seed_r0(:, n), &
                  init_conds%seed_r1(:, n), &
                  init_conds%seed_density(n), init_conds%seed_density2(n), NDIM, &
                  init_conds%seed_width(n), &
                  init_conds%seed_falloff(n))

             ! Add electrons and/or ions depending on the seed charge type
             ! (positive, negative or neutral)
             if (init_conds%seed_charge_type(n) <= 0) then
                box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + density
             end if

             if (init_conds%seed_charge_type(n) >= 0) then
                box%cc(IJK, i_1pos_ion) = box%cc(IJK, i_1pos_ion) + density
             end if
          end do
       end do
    end do
  end subroutine my_init_cond

end module m_user
