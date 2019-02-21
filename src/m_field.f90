#include "../afivo/src/cpp_macros.h"
!> Module to compute electric fields
module m_field
  use m_af_all
  use m_types
  use m_streamer

  implicit none
  private

  !> Use a table with fields versus time
  logical :: field_table_use

  !> List of times
  real(dp), allocatable :: field_table_times(:)

  !> List of fields
  real(dp), allocatable :: field_table_fields(:)

  !> Start modifying the vertical background field after this time
  real(dp) :: field_mod_t0 = 1e99_dp

  !> Stop modifying the vertical background field after this time
  real(dp) :: field_mod_t1 = 1e99_dp

  !> Add this uniform background field (V/m)
  real(dp) :: field_background(NDIM) = 0.0_dp

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
  real(dp), public, protected :: field_voltage

  logical  :: field_stability_search    = .false.
  real(dp) :: field_stability_zmin      = 0.2_dp
  real(dp) :: field_stability_zmax      = 1.0_dp
  real(dp) :: field_stability_threshold = 3e6_dp

  real(dp) :: field_point_charge = 0.0_dp
#if NDIM == 2
  real(dp) :: field_point_r0(NDIM) = [0.0_dp, -1.0_dp]
#elif NDIM == 3
  real(dp) :: field_point_r0(NDIM) = [0.0_dp, 0.0_dp, -1.0_dp]
#endif

  character(string_len) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_set_rhs
  public :: field_from_potential
  public :: field_get_amplitude
  public :: field_set_voltage

  public :: field_bc_homogeneous

contains

  !> Initialize this module
  subroutine field_initialize(tree, cfg, mg)
    use m_config
    use m_table_data
    use m_user_methods
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg !< Settings
    type(mg_t), intent(inout)  :: mg  !< Multigrid option struct
    character(len=string_len)  :: field_table

    field_table = undefined_str
    call CFG_add_get(cfg, "field_table", field_table, &
         "File containing applied electric field (V/m) versus time")
    if (field_table /= undefined_str) then
       field_table_use = .true.
       call table_from_file(field_table, "field_vs_time", &
            field_table_times, field_table_fields)
    else
       field_table_use = .false.
    end if

    call CFG_add_get(cfg, "field_mod_t0", field_mod_t0, &
         "Modify electric field after this time (s)")
    call CFG_add_get(cfg, "field_mod_t1", field_mod_t1, &
         "Modify electric field up to this time (s)")
    call CFG_add_get(cfg, "field_background", field_background, &
         "Add this uniform background field (V/m)")
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

    call CFG_add_get(cfg, "field_stability_search", field_stability_search, &
         "If true, enable mode to search stability field")
    call CFG_add_get(cfg, "field_stability_zmin", field_stability_zmin, &
         "Start lowering background field above this relative position")
    call CFG_add_get(cfg, "field_stability_zmax", field_stability_zmax, &
         "At this relative position the background field will be zero")
    call CFG_add_get(cfg, "field_stability_threshold", field_stability_threshold, &
         "Use location of maximal field if above this threshold (V/m)")

    call CFG_add_get(cfg, "field_point_charge", field_point_charge, &
         "Charge (in C) of point charge")
    call CFG_add_get(cfg, "field_point_r0", field_point_r0, &
         "Relative position of point charge (outside domain)")

    field_voltage = -ST_domain_len(NDIM) * field_amplitude

    if (associated(user_potential_bc)) then
       mg%sides_bc => user_potential_bc
    else
       ! Use one of the predefined boundary conditions
       select case (field_bc_type)
       case ("homogeneous")
          mg%sides_bc => field_bc_homogeneous
       case ("point_charge")
          mg%sides_bc => field_bc_point_charge
       case default
          error stop "field_bc_select error: invalid condition"
       end select
    end if

    ! Set the multigrid options. First define the variables to use
    mg%i_phi = i_phi
    mg%i_tmp = i_tmp
    mg%i_rhs = i_rhs
    if (ST_use_dielectric) mg%i_eps = i_eps

    ! This automatically handles cylindrical symmetry
    mg%box_op => mg_auto_op
    mg%box_gsrb => mg_auto_gsrb
    mg%box_corr => mg_auto_corr
    mg%box_stencil => mg_auto_stencil
    mg%sides_rb => mg_sides_rb

    call af_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)
    call af_set_cc_methods(tree, i_electric_fld, &
         af_bc_neumann_zero, af_gc_interp)

  end subroutine field_initialize

  subroutine field_set_rhs(tree, s_in)
    use m_units_constants
    use m_chemistry
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: s_in
    real(dp), parameter       :: fac = -UC_elem_charge / UC_eps0
    real(dp)                  :: q
    integer                   :: lvl, i, id, nc, n, ix

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id, n, ix, q)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)

          tree%boxes(id)%cc(DTIMES(:), i_rhs) = 0.0_dp

          do n = 1, size(charged_species_itree)
             ix = charged_species_itree(n) + s_in
             q = charged_species_charge(n) * fac

             tree%boxes(id)%cc(DTIMES(:), i_rhs) = &
                  tree%boxes(id)%cc(DTIMES(:), i_rhs) + &
                  q * tree%boxes(id)%cc(DTIMES(:), ix)
          end do
       end do
       !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine field_set_rhs

  !> Compute electric field on the tree. First perform multigrid to get electric
  !> potential, then take numerical gradient to geld field.
  subroutine field_compute(tree, mg, s_in, time, have_guess)
    use m_units_constants
    use m_chemistry
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg ! Multigrid option struct
    integer, intent(in)       :: s_in
    real(dp), intent(in)      :: time
    logical, intent(in)       :: have_guess
    integer                   :: i

    call field_set_rhs(tree, s_in)
    call field_set_voltage(time)

    if (.not. have_guess) then
       ! Perform a FMG cycle when we have no guess
       call mg_fas_fmg(tree, mg, .false., have_guess)
    else
       ! Perform cheaper V-cycles
       do i = 1, ST_multigrid_num_vcycles
          call mg_fas_vcycle(tree, mg, .false.)
       end do
    end if

    ! Compute field from potential
    call af_loop_box(tree, field_from_potential)

    ! Set the field norm also in ghost cells
    call af_gc_tree(tree, i_electric_fld)
  end subroutine field_compute

  !> Compute the electric field at a given time
  function field_get_amplitude(time) result(electric_fld)
    use m_units_constants
    use m_lookup_table
    real(dp), intent(in)   :: time
    real(dp)               :: electric_fld, t_rel
    type(af_loc_t)         :: loc_field
    real(dp)               :: r(NDIM), zrel, max_fld

    if (field_table_use) then
       call LT_lin_interp_list(field_table_times, field_table_fields, &
            time, electric_fld)
    else
       ! TODO: simplify stuff below
       t_rel = time - field_mod_t0
       t_rel = min(t_rel, field_mod_t1-field_mod_t0)

       if (t_rel > 0) then
          electric_fld = field_amplitude * exp(-t_rel/field_decay_time) + &
               t_rel * field_lin_deriv + &
               field_sin_amplitude * &
               sin(t_rel * field_sin_freq * 2 * UC_pi)
       else
          electric_fld = field_amplitude
       end if
    end if

  end function field_get_amplitude

  !> Compute the voltage at a given time
  subroutine field_set_voltage(time)
    real(dp), intent(in) :: time

    field_voltage = -ST_domain_len(NDIM) * field_get_amplitude(time)
  end subroutine field_set_voltage

  !> This fills ghost cells near physical boundaries for the potential
  subroutine field_bc_homogeneous(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (af_neighb_dim(nb) == NDIM) then
       if (af_neighb_low(nb)) then
          bc_type = af_bc_dirichlet
          bc_val = 0.0_dp
       else
          bc_type = af_bc_dirichlet
          bc_val  = field_voltage
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine field_bc_homogeneous

  !> Create a field of the form E = E_0 - c / r^2
  subroutine field_bc_point_charge(box, nb, iv, coords, bc_val, bc_type)
    use m_units_constants
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: n
    real(dp)                :: r0(NDIM), r1(NDIM), q

    bc_type  = af_bc_dirichlet
    q        = field_point_charge / (4 * UC_pi * UC_eps0)
    r0       = field_point_r0 * ST_domain_len
    r1       = r0
    r1(NDIM) = -r0(NDIM)

    if (ST_cylindrical .and. nb == af_neighb_lowx) then
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    else
       bc_type = af_bc_dirichlet
       do n = 1, box%n_cell**(NDIM-1)
          bc_val(n) = q/norm2(coords(:, n) - r0) - &
               q/norm2(coords(:, n) - r1)
       end do
    end if
  end subroutine field_bc_point_charge

  !> Compute electric field from electrical potential
  subroutine field_from_potential(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    real(dp)                   :: inv_dr(NDIM)

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if NDIM == 2
    box%fc(1:nc+1, 1:nc, 1, electric_fld) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi)) + &
         field_background(1)
    box%fc(1:nc, 1:nc+1, 2, electric_fld) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi)) + &
         field_background(2)

    if (ST_use_dielectric) then
       ! Compute fields at the boundaries of the box, where eps can change
       box%fc(1, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
            (box%cc(0, 1:nc, i_phi) - box%cc(1, 1:nc, i_phi)) * &
            box%cc(0, 1:nc, i_eps) / &
            (box%cc(1, 1:nc, i_eps) + box%cc(0, 1:nc, i_eps))
       box%fc(nc+1, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
            (box%cc(nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, i_phi)) * &
            box%cc(nc+1, 1:nc, i_eps) / &
            (box%cc(nc+1, 1:nc, i_eps) + box%cc(nc, 1:nc, i_eps))
       box%fc(1:nc, 1, 2, electric_fld) = 2 * inv_dr(2) * &
            (box%cc(1:nc, 0, i_phi) - box%cc(1:nc, 1, i_phi)) * &
            box%cc(1:nc, 0, i_eps) / &
            (box%cc(1:nc, 1, i_eps) + box%cc(1:nc, 0, i_eps))
       box%fc(1:nc, nc+1, 2, electric_fld) = 2 * inv_dr(2) * &
            (box%cc(1:nc, nc, i_phi) - box%cc(1:nc, nc+1, i_phi)) * &
            box%cc(1:nc, nc+1, i_eps) / &
            (box%cc(1:nc, nc+1, i_eps) + box%cc(1:nc, nc, i_eps))
    end if

    box%cc(1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 2, electric_fld))**2)
#elif NDIM == 3
    box%fc(1:nc+1, 1:nc, 1:nc, 1, electric_fld) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, 1:nc, i_phi) - &
         box%cc(1:nc+1, 1:nc, 1:nc, i_phi)) + &
         field_background(1)
    box%fc(1:nc, 1:nc+1, 1:nc, 2, electric_fld) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, 1:nc, i_phi) - &
         box%cc(1:nc, 1:nc+1, 1:nc, i_phi)) + &
         field_background(2)
    box%fc(1:nc, 1:nc, 1:nc+1, 3, electric_fld) = inv_dr(3) * &
         (box%cc(1:nc, 1:nc, 0:nc, i_phi) - &
         box%cc(1:nc, 1:nc, 1:nc+1, i_phi)) + &
         field_background(3)

    if (ST_use_dielectric) then
       ! Compute fields at the boundaries of the box, where eps can change
       box%fc(1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
            (box%cc(0, 1:nc, 1:nc, i_phi) - box%cc(1, 1:nc, 1:nc, i_phi)) * &
            box%cc(0, 1:nc, 1:nc, i_eps) / &
            (box%cc(1, 1:nc, 1:nc, i_eps) + box%cc(0, 1:nc, 1:nc, i_eps))
       box%fc(nc+1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
            (box%cc(nc, 1:nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, 1:nc, i_phi)) * &
            box%cc(nc+1, 1:nc, 1:nc, i_eps) / &
            (box%cc(nc+1, 1:nc, 1:nc, i_eps) + box%cc(nc, 1:nc, 1:nc, i_eps))
       box%fc(1:nc, 1, 1:nc, 2, electric_fld) = 2 * inv_dr(2) * &
            (box%cc(1:nc, 0, 1:nc, i_phi) - box%cc(1:nc, 1, 1:nc, i_phi)) * &
            box%cc(1:nc, 0, 1:nc, i_eps) / &
            (box%cc(1:nc, 1, 1:nc, i_eps) + box%cc(1:nc, 0, 1:nc, i_eps))
       box%fc(1:nc, nc+1, 1:nc, 2, electric_fld) = 2 * inv_dr(2) * &
            (box%cc(1:nc, nc, 1:nc, i_phi) - box%cc(1:nc, nc+1, 1:nc, i_phi)) * &
            box%cc(1:nc, nc+1, 1:nc, i_eps) / &
            (box%cc(1:nc, nc+1, 1:nc, i_eps) + box%cc(1:nc, nc, 1:nc, i_eps))
       box%fc(1:nc, 1:nc, 1, 3, electric_fld) = 2 * inv_dr(3) * &
            (box%cc(1:nc, 1:nc, 0, i_phi) - box%cc(1:nc, 1:nc, 1, i_phi)) * &
            box%cc(1:nc, 1:nc, 0, i_eps) / &
            (box%cc(1:nc, 1:nc, 1, i_eps) + box%cc(1:nc, 1:nc, 0, i_eps))
       box%fc(1:nc, 1:nc, nc+1, 3, electric_fld) = 2 * inv_dr(3) * &
            (box%cc(1:nc, 1:nc, nc, i_phi) - box%cc(1:nc, 1:nc, nc+1, i_phi)) * &
            box%cc(1:nc, 1:nc, nc+1, i_eps) / &
            (box%cc(1:nc, 1:nc, nc+1, i_eps) + box%cc(1:nc, 1:nc, nc, i_eps))
    end if

    box%cc(1:nc, 1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 1:nc, 2, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 1:nc, 3, electric_fld) + &
         box%fc(1:nc, 1:nc, 2:nc+1, 3, electric_fld))**2)
#endif

  end subroutine field_from_potential

end module m_field
