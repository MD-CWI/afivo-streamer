module m_field_$Dd
  use m_streamer
  use m_a$D_all

  implicit none
  private

  ! Start modifying the vertical background field after this time
  real(dp), protected :: ST_electric_fld_y_mod_t0

  ! Amplitude of sinusoidal modification
  real(dp), protected :: ST_electric_fld_y_sin_amplitude

  ! Frequency (Hz) of sinusoidal modification
  real(dp), protected :: ST_electric_fld_y_sin_freq

  ! Linear derivative of background field
  real(dp), protected :: ST_electric_fld_y_lin_deriv

  ! Decay time of background field
  real(dp), protected :: ST_electric_fld_y_decay

  ! Start modifying the horizontal background field after this time
  real(dp), protected :: ST_electric_fld_x_mod_t0

  ! Amplitude of sinusoidal modification
  real(dp), protected :: ST_electric_fld_x_sin_amplitude

  ! Frequency (Hz) of sinusoidal modification
  real(dp), protected :: ST_electric_fld_x_sin_freq

  ! Linear derivative of background field
  real(dp), protected :: ST_electric_fld_x_lin_deriv

  ! Decay time of background field
  real(dp), protected :: ST_electric_fld_x_decay

    ! The applied electric field (vertical direction)
  real(dp), protected :: ST_applied_electric_fld_y

  ! The applied electric field (horizontal direction)
  real(dp), protected :: ST_applied_electric_fld_x

  ! The applied voltage (vertical direction)
  real(dp), protected :: ST_applied_voltage

  public :: field_bc_select
  public :: field_from_potential

contains

  !> Set boundary conditions for the electric potential/field
  subroutine field_initialize(cfg, mg)
    type(CFG_t), intent(inout)  :: cfg !< Settings
    type(mg$D_t), intent(inout) :: mg  !< Multigrid option struct

    call CFG_add(ST_config, "electric_fld_y_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time (s)")
    call CFG_add(ST_config, "electric_fld_y_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add(ST_config, "electric_fld_y_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add(ST_config, "electric_fld_y_lin_deriv", 0.0_dp, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add(ST_config, "electric_fld_y_decay", huge(1.0_dp), &
         "Decay time of field (s)")

    call CFG_add(cfg, "applied_electric_fld_y", applied_electric_fld_y.0d7, &
         "The applied electric field (V/m) (vertical)")
    call CFG_add(cfg, "applied_electric_fld_x", 0.0d0, &
         "The applied electric field (V/m) (horizontal)")

    call CFG_add(ST_config, "electric_fld_x_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time (s)")
    call CFG_add(ST_config, "electric_fld_x_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add(ST_config, "electric_fld_x_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add(ST_config, "electric_fld_x_lin_deriv", 0.0_dp, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add(ST_config, "electric_fld_x_decay", huge(1.0_dp), &
         "Decay time of field (s)")

    ST_applied_voltage = -ST_domain_len * ST_applied_electric_fld_y

    select case (ST_field_bc)
    case ("homogeneous")
       mg%sides_bc => field_bc_homogeneous
    case default
       error stop "field_bc_select error: invalid condition"
    end select
  end subroutine field_initialize

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
       box%cc(1:nc, nc+1, iv) = ST_applied_voltage
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
       box%cc(1:nc, 1:nc, nc+1, iv) = ST_applied_voltage
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
