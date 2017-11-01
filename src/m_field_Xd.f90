#include "../afivo/src/cpp_macros_$Dd.h"
module m_field_$Dd
  use m_streamer
  use m_a$D_all

  implicit none
  private

  !> Start modifying the vertical background field after this time
  real(dp) :: field_mod_t0 = 1e99_dp

  !> Amplitude of sinusoidal modification
  real(dp) :: field_sin_amplitude = 0.0_dp

  !> Frequency (Hz) of sinusoidal modification
  real(dp) :: field_sin_freq = 0.0_dp

  !> Linear derivative of background field
  real(dp) :: field_lin_deriv = 0.0_dp

  !> Decay time of background field
  real(dp) :: field_decay_time = huge(1.0_dp)

  !> The applied electric field (vertical direction)
  real(dp) :: field_amplitude = 1.0e6_dp

  !> The applied voltage (vertical direction)
  real(dp) :: field_voltage

  !> Drop-off radius
  real(dp) :: field_dropoff_radius = 1e-3_dp

  !> Relative width over which the potential drops
  real(dp) :: field_dropoff_relwidth = 0.5_dp

  !> Location from which the field drops off (set below)
  real(dp) :: field_dropoff_pos(2) = 0.0_dp

  character(ST_slen) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_from_potential
  public :: field_get_amplitude
  public :: field_set_voltage

  public :: field_bc_homogeneous
  public :: field_bc_dropoff_lin
  public :: field_bc_dropoff_log

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

    call CFG_add_get(cfg, "field_dropoff_radius", field_dropoff_radius, &
         "Potential stays constant up to this radius")
    call CFG_add_get(cfg, "field_dropoff_relwidth", field_dropoff_relwidth, &
         "Relative width over which the potential drops")

    field_voltage = -ST_domain_len * field_amplitude

    select case (field_bc_type)
    case ("homogeneous")
       mg%sides_bc => field_bc_homogeneous
    case ("dropoff_lin")
       if (ST_cylindrical) then
          field_dropoff_pos(:) = 0.0_dp
       else
          field_dropoff_pos(:) = 0.5_dp
       end if

       mg%sides_bc => field_bc_dropoff_lin
    case ("dropoff_log")
       if (ST_cylindrical) then
          field_dropoff_pos(:) = 0.0_dp
       else
          field_dropoff_pos(:) = 0.5_dp
       end if

       mg%sides_bc => field_bc_dropoff_log
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

    if (.not. have_guess) then
       call mg$D_fas_fmg(tree, mg, .false., have_guess)
       call mg$D_fas_fmg(tree, mg, .false., .true.)
    else
       call mg$D_fas_vcycle(tree, mg, .false.)
       call mg$D_fas_vcycle(tree, mg, .false.)
    end if

    ! Compute field from potential
    call a$D_loop_box(tree, field_from_potential)

    ! Set the fields also in ghost cells
    call a$D_gc_tree(tree, i_electric_fld, a$D_gc_interp, a$D_bc_neumann_zero)
    call a$D_gc_tree(tree, i_Ex, a$D_gc_interp, a$D_bc_neumann_zero)
    call a$D_gc_tree(tree, i_Ey, a$D_gc_interp, a$D_bc_neumann_zero)
#if $D == 3
    call a$D_gc_tree(tree, i_Ez, a$D_gc_interp, a$D_bc_neumann_zero)
#endif
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
    use m_a$D_ghostcell
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
#if $D == 2
    case (a$D_neighb_lowy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc,    0, iv) = 0
    case (a$D_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = field_voltage
#elif $D == 3
    case (a3_neighb_lowz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc,    0, iv) = 0
    case (a3_neighb_highz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc, nc+1, iv) = field_voltage
#endif
    case default
       call a$D_bc_continuous(box, nb, iv, bc_type)
    end select

  end subroutine field_bc_homogeneous

  subroutine field_bc_dropoff_lin(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc, i
#if $D == 3
    integer                      :: j
#endif
    real(dp)                     :: rr($D), rdist

    nc = box%n_cell

    select case (nb)
#if $D == 2
    case (a2_neighb_highy)
       bc_type = af_bc_dirichlet

       do i = 1, nc
          rr = a2_r_cc(box, [i, 0])
          rdist = abs(rr(1) - field_dropoff_pos(1))
          rdist = (rdist - field_dropoff_radius) / &
               (field_dropoff_relwidth * ST_domain_len)

          if (rdist < 0) then
             box%cc(i, nc+1, iv) = field_voltage
          else
             box%cc(i, nc+1, iv) = field_voltage * &
                  max(0.0_dp, (1 - rdist))
          end if
       end do
#elif $D == 3
    case (a3_neighb_highz)
       bc_type = af_bc_dirichlet

       do j = 1, nc
          do i = 1, nc
             rr = a3_r_cc(box, [i, j, 0])
             rdist = norm2(rr(1:2) - field_dropoff_pos(1:2))
             rdist = (rdist - field_dropoff_radius) / &
                  (field_dropoff_relwidth * ST_domain_len)

             if (rdist < 0) then
                box%cc(i, j, nc+1, iv) = field_voltage
             else
                box%cc(i, j, nc+1, iv) = field_voltage * &
                     max(0.0_dp, (1 - rdist))
             end if
          end do
       end do
#endif
    case default
       call field_bc_homogeneous(box, nb, iv, bc_type)
    end select
  end subroutine field_bc_dropoff_lin

  subroutine field_bc_dropoff_log(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc, i
#if $D == 3
    integer                      :: j
#endif
    real(dp)                     :: rr($D), rdist, tmp

    nc = box%n_cell
    tmp = field_dropoff_relwidth * ST_domain_len

    select case (nb)
#if $D == 2
    case (a2_neighb_highy)
       bc_type = af_bc_dirichlet

       do i = 1, nc
          rr = a2_r_cc(box, [i, 0])
          rdist = abs(rr(1) - field_dropoff_pos(1))

          if (rdist < field_dropoff_radius) then
             box%cc(i, nc+1, iv) = field_voltage
          else
             box%cc(i, nc+1, iv) = field_voltage * &
                  log(1 + tmp/rdist) / log(1 + tmp/field_dropoff_radius)
          end if
       end do
#elif $D == 3
    case (a3_neighb_highz)
       bc_type = af_bc_dirichlet

       do j = 1, nc
          do i = 1, nc
             rr = a3_r_cc(box, [i, j, 0])
             rdist = norm2(rr(1:2) - field_dropoff_pos(1:2))

             if (rdist < field_dropoff_radius) then
                box%cc(i, j, nc+1, iv) = field_voltage
             else
                box%cc(i, j, nc+1, iv) = field_voltage * &
                     log(1 + tmp/rdist) / log(1 + tmp/field_dropoff_radius)
             end if
          end do
       end do
#endif
    case default
       call field_bc_homogeneous(box, nb, iv, bc_type)
    end select
  end subroutine field_bc_dropoff_log

  !> Compute electric field from electrical potential
  subroutine field_from_potential(box)
    type(box$D_t), intent(inout) :: box
    integer                     :: nc, IJK
    real(dp)                    :: fac

    nc     = box%n_cell
    fac = 0.5_dp / box%dr

    do KJI_DO(1,nc)
#if $D == 2
       box%cc(i, j, i_Ex) = (box%cc(i-1, j, i_phi) - &
            box%cc(i+1, j, i_phi)) * fac
       box%cc(i, j, i_Ey) = (box%cc(i, j-1, i_phi) - &
            box%cc(i, j+1, i_phi)) * fac
       box%cc(i, j, i_electric_fld) = norm2(box%cc(i, j, i_Ex:i_Ey))
#elif $D == 3
       box%cc(i, j, k, i_Ex) = (box%cc(i-1, j, k, i_phi) - &
            box%cc(i+1, j, k, i_phi)) * fac
       box%cc(i, j, k, i_Ey) = (box%cc(i, j-1, k, i_phi) - &
            box%cc(i, j+1, k, i_phi)) * fac
       box%cc(i, j, k, i_Ez) = (box%cc(i, j, k-1, i_phi) - &
            box%cc(i, j, k+1, i_phi)) * fac
       box%cc(i, j, k, i_electric_fld) = &
            norm2(box%cc(i, j, k, i_Ex:i_Ez))
#endif
    end do; CLOSE_DO
  end subroutine field_from_potential

end module m_field_$Dd
