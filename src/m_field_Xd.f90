#include "../afivo/src/cpp_macros_$Dd.h"
module m_field_$Dd
  use m_streamer
  use m_a$D_all

  implicit none
  private

  ! Start modifying the vertical background field after this time
  real(dp), protected :: field_mod_t0 = 1e99_dp

  ! Amplitude of sinusoidal modification
  real(dp), protected :: field_sin_amplitude = 0.0_dp

  ! Frequency (Hz) of sinusoidal modification
  real(dp), protected :: field_sin_freq = 0.0_dp

  ! Linear derivative of background field
  real(dp), protected :: field_lin_deriv = 0.0_dp

  ! Decay time of background field
  real(dp), protected :: field_decay_time = huge(1.0_dp)

  ! The applied electric field (vertical direction)
  real(dp), protected :: field_amplitude = 1.0e6_dp

  ! The applied voltage (vertical direction)
  real(dp), protected :: field_voltage

  character(ST_slen) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_from_potential
  public :: field_get_amplitude
  public :: field_set_voltage

  public :: field_bc_homogeneous

contains

  !> Set boundary conditions for the electric potential/field
  subroutine field_initialize(cfg, mg)
    type(CFG_t), intent(inout)  :: cfg !< Settings
    type(mg$D_t), intent(inout) :: mg  !< Multigrid option struct

    call CFG_add_get(cfg, "field_mod_t0", field_mod_t0, &
         "Modify electric field after this time (s)")
    call CFG_add_get(cfg, "field_sin_amplitude", field_sin_amplitude, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add_get(cfg, "field_sin_freq", field_sin_freq, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add_get(cfg, "field_lin_deriv", field_lin_deriv, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add_get(cfg, "field_decay_time", field_decay_time, &
         "Decay time of field (s)")
    call CFG_add_get(cfg, "field_amplitude", field_amplitude, &
         "The applied electric field (V/m) (vertical)")
    call CFG_add_get(cfg, "field_bc_type", field_bc_type, &
         "Type of boundary condition to use (homogeneous, ...)")

    field_voltage = -ST_domain_len * field_amplitude

    select case (field_bc_type)
    case ("homogeneous")
       mg%sides_bc => field_bc_homogeneous
    case default
       error stop "field_bc_select error: invalid condition"
    end select
  end subroutine field_initialize

  !> Compute electric field on the tree. First perform multigrid to get electric
  !> potential, then take numerical gradient to geld field.
  subroutine field_compute(tree, mg, have_guess)
    use m_units_constants
    type(a$D_t), intent(inout) :: tree
    type(mg$D_t), intent(in)   :: mg ! Multigrid option struct
    logical, intent(in)        :: have_guess
    real(dp), parameter        :: fac = UC_elem_charge / UC_eps0
    integer                    :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(DTIMES(:), i_rhs) = fac * (&
               tree%boxes(id)%cc(DTIMES(:), i_electron) - &
               tree%boxes(id)%cc(DTIMES(:), i_pos_ion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    call field_set_voltage(ST_time)
    call mg$D_fas_fmg(tree, mg, .false., have_guess)

    ! Compute field from potential
    call a$D_loop_box(tree, field_from_potential)

    ! Set the field norm also in ghost cells
    call a$D_gc_tree(tree, i_electric_fld, a$D_gc_interp, a$D_bc_neumann_zero)
  end subroutine field_compute

  !> Compute the electric field at a given time
  function field_get_amplitude(time) result(electric_fld)
    use m_units_constants
    real(dp), intent(in) :: time
    real(dp)             :: electric_fld, t_rel

    t_rel = time - field_mod_t0
    if (t_rel > 0) then
       electric_fld = field_amplitude * exp(-t_rel/field_decay_time) + &
            t_rel * field_lin_deriv + &
            field_sin_amplitude * &
            sin(t_rel * field_sin_freq * 2 * UC_pi)
    else
       electric_fld = field_amplitude
    end if
  end function field_get_amplitude

  !> Compute the voltage at a given time
  subroutine field_set_voltage(time)
    real(dp), intent(in) :: time

    field_voltage = -ST_domain_len * field_get_amplitude(time)
  end subroutine field_set_voltage

  !> This fills ghost cells near physical boundaries for the potential
  subroutine field_bc_homogeneous(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
#if $D == 2
    case (a$D_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, iv) = 0
    case (a$D_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0
    case (a$D_neighb_lowy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc,    0, iv) = 0
    case (a$D_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = field_voltage
#elif $D == 3
    case (a3_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, 1:nc, iv) = 0
    case (a3_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, 1:nc, iv) = 0
    case (a3_neighb_lowy)
       bc_type = af_bc_neumann
       box%cc(1:nc,    0, 1:nc, iv) = 0
    case (a3_neighb_highy)
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, 1:nc, iv) = 0
    case (a3_neighb_lowz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc,    0, iv) = 0
    case (a3_neighb_highz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc, nc+1, iv) = field_voltage
#endif
    end select

  end subroutine field_bc_homogeneous

  !> Compute electric field from electrical potential
  subroutine field_from_potential(box)
    type(box$D_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if $D == 2
    box%fc(1:nc+1, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 2, electric_fld))**2)
#elif $D == 3
    box%fc(1:nc+1, 1:nc, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, 1:nc, i_phi) - &
         box%cc(1:nc+1, 1:nc, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 1:nc, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, 1:nc, i_phi) - &
         box%cc(1:nc, 1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc, 1:nc+1, 3, electric_fld) = inv_dr * &
         (box%cc(1:nc, 1:nc, 0:nc, i_phi) - &
         box%cc(1:nc, 1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 1:nc, 2, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 1:nc, 3, electric_fld) + &
         box%fc(1:nc, 1:nc, 2:nc+1, 3, electric_fld))**2)
#endif

  end subroutine field_from_potential

end module m_field_$Dd
