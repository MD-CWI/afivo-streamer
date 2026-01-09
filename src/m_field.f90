!> Module to compute electric fields
module m_field
#include "../afivo/src/cpp_macros.h"
  use m_af_all
  use m_types
  use m_streamer

  implicit none
  private

  integer, parameter :: scalar_voltage = 1
  integer, parameter :: tabulated_voltage = 2

  !> How the electric field or voltage is specified
  integer :: field_given_by = -1

  !> List of times
  real(dp), allocatable :: field_table_times(:)

  !> List of voltages
  real(dp), allocatable :: field_table_values(:)

  !> Linear rise time of field (s)
  real(dp) :: field_rise_time = 0.0_dp

  !> Pulse width excluding rise and fall time
  real(dp) :: field_pulse_width = huge_real

  !> Number of voltage pulses
  integer :: field_num_pulses = 1

  !> Time of one complete voltage pulse
  real(dp), public, protected :: field_pulse_period = huge_real

  !> The (initial) vertical applied electric field
  real(dp) :: field_amplitude = undefined_real

  !> The applied voltage (vertical direction)
  real(dp) :: field_voltage = undefined_real

  !> Whether electrode 1 is grounded or at the applied voltage
  logical :: field_electrode_grounded = .false.

  !> Whether electrode 2 is grounded or at the applied voltage
  logical :: field_electrode2_grounded = .false.

  !> Electrode 1: first relative coordinate
  real(dp) :: rod_r0(NDIM) = -1.0e100_dp

  !> Electrode 1: second relative coordinate
  real(dp) :: rod_r1(NDIM) = -1.0e100_dp

  !> Electrode 2: first relative coordinate
  real(dp) :: rod2_r0(NDIM) = -1.0e100_dp

  !> Electrode 2: second relative coordinate
  real(dp) :: rod2_r1(NDIM) = -1.0e100_dp

  !> Electrode 1 radius (in m)
  real(dp) :: rod_radius = -1.0e100_dp

  !> Electrode 2 radius (in m)
  real(dp) :: rod2_radius = -1.0e100_dp

  !> Electrode 1: fraction of conical part (if conical)
  real(dp) :: cone_length_frac = -1.0e100_dp

  !> Electrode 2: fraction of conical part (if conical)
  real(dp) :: cone2_length_frac = -1.0e100_dp

  !> Electrode 1: tip radius (if conical)
  real(dp) :: cone_tip_radius = -1.0e100_dp

  !> Electrode 2: tip radius (if conical)
  real(dp) :: cone2_tip_radius = -1.0e100_dp

  ! Internal variables

  !> Electrode 1: 'origin' of spherical tip (if conical)
  real(dp) :: cone_tip_center(NDIM)
  !> Electrode 2: 'origin' of spherical tip (if conical)
  real(dp) :: cone2_tip_center(NDIM)

  !> Electrode 1: radius of curvature of spherical tip (if conical)
  real(dp) :: cone_tip_r_curvature
  !> Electrode 2: radius of curvature of spherical tip (if conical)
  real(dp) :: cone2_tip_r_curvature

  !> The current applied voltage
  real(dp), public, protected :: current_voltage = 0.0_dp

  character(string_len) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_set_rhs
  public :: field_set_voltage

  public :: field_bc_homogeneous
  public :: field_from_potential
  public :: field_compute_energy
  public :: field_get_E_vector

contains

  !> Initialize this module
  subroutine field_initialize(tree, cfg, mg)
    use m_config
    use m_table_data
    use m_user_methods
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg !< Settings
    type(mg_t), intent(inout)  :: mg  !< Multigrid option struct
    character(len=string_len)  :: given_by, user_value
    character(len=string_len)  :: electrode_type
    integer                    :: first_blank

    ! This is for backward compatibility
    call CFG_add_get(cfg, "field_amplitude", field_amplitude, &
         "The (initial) vertical applied electric field (V/m)")

    given_by = undefined_str
    call CFG_add_get(cfg, "field_given_by", given_by, &
         "How the electric field or voltage is specified")

    ! Split given_by string
    given_by = adjustl(given_by)
    first_blank = index(given_by, " ")
    user_value = adjustl(given_by(first_blank:))
    given_by = given_by(1:first_blank-1)

    select case (given_by)
    case ("voltage")
       field_given_by = scalar_voltage
       read(user_value, *) field_voltage
    case ("field")
       field_given_by = scalar_voltage
       read(user_value, *) field_voltage
       ! Convert to a voltage
       field_voltage = -ST_domain_len(NDIM) * field_voltage
    case ("voltage_table")
       field_given_by = tabulated_voltage

       ! Read in the table
       call table_from_file(trim(user_value), "voltage_vs_time", &
            field_table_times, field_table_values)
    case ("field_table")
       field_given_by = tabulated_voltage

       ! Read in the table
       call table_from_file(trim(user_value), "field_vs_time", &
            field_table_times, field_table_values)
       ! Convert to a voltage
       field_table_values = -ST_domain_len(NDIM) * field_table_values
    case (undefined_str)
       if (field_amplitude <= undefined_real) then
          error stop "field_amplitude not specified"
       else
          print *, "Warning: field_amplitude is deprecated, use field_given_by"
          field_given_by = scalar_voltage
          field_voltage = -ST_domain_len(NDIM) * field_amplitude
       end if
    case default
       print *, "field_given_by value: ", trim(given_by), ", options are:"
       print *, "1. voltage <value in V>"
       print *, "2. field <value in V/m>"
       print *, "3. voltage_table <filename>"
       print *, "4. field_table <filename>"
       error stop "Unknown field_given_by value"
    end select

    call CFG_add_get(cfg, "field_rise_time", field_rise_time, &
         "Linear rise time of field (s)")
    call CFG_add_get(cfg, "field_pulse_width", field_pulse_width, &
         "Pulse width excluding rise and fall time (s)")
    call CFG_add_get(cfg, "field_num_pulses", field_num_pulses, &
         "Number of voltage pulses (default: 1)")
    call CFG_add_get(cfg, "field_pulse_period", field_pulse_period, &
         "Time of one complete voltage pulse (s)")

    if (field_pulse_width < huge_real .and. field_rise_time <= 0) &
         error stop "Set field_rise_time when using field_pulse_width"

    if (field_num_pulses > 1) then
       if (field_pulse_period >= huge_real) &
            error stop "field_num_pulses > 1 requires field_pulse_period"
       if (field_pulse_width >= huge_real) &
            error stop "field_num_pulses > 1 requires field_pulse_width"
       if (field_pulse_width + 2 * field_rise_time > field_pulse_period) &
            error stop "field_pulse_period shorter than one full pulse"
    end if

    call CFG_add_get(cfg, "field_bc_type", field_bc_type, &
         "Boundary condition for electric potential")

    !< [electrode_settings]
    call CFG_add_get(cfg, "field_electrode_grounded", field_electrode_grounded, &
         "Whether electrode 1 is grounded or at the applied voltage")
    call CFG_add_get(cfg, "field_electrode2_grounded", field_electrode2_grounded, &
         "Whether electrode 2 is grounded or at the applied voltage")
    call CFG_add_get(cfg, "field_rod_r0", rod_r0, &
         "Electrode 1: first relative coordinate")
    call CFG_add_get(cfg, "field_rod_r1", rod_r1, &
         "Electrode 1: second relative coordinate")
    call CFG_add_get(cfg, "field_rod2_r0", rod2_r0, &
         "Electrode 2: first relative coordinate")
    call CFG_add_get(cfg, "field_rod2_r1", rod2_r1, &
         "Electrode 2: second relative coordinate")
    call CFG_add_get(cfg, "field_rod_radius", rod_radius, &
         "Electrode 1 radius (in m)")
    call CFG_add_get(cfg, "field_rod2_radius", rod2_radius, &
         "Electrode 2 radius (in m)")
    call CFG_add_get(cfg, "cone_tip_radius", cone_tip_radius, &
         "Electrode 1: tip radius (if conical)")
    call CFG_add_get(cfg, "cone_length_frac", cone_length_frac, &
         "Electrode 1: fraction of conical part (if conical)")
    call CFG_add_get(cfg, "cone2_tip_radius", cone2_tip_radius, &
         "Electrode 2: tip radius (if conical)")
    call CFG_add_get(cfg, "cone2_length_frac", cone2_length_frac, &
         "Electrode 2: fraction of conical part (if conical)")

    rod_r0 = ST_domain_origin + rod_r0 * ST_domain_len
    rod_r1 = ST_domain_origin + rod_r1 * ST_domain_len
    rod2_r0 = ST_domain_origin + rod2_r0 * ST_domain_len
    rod2_r1 = ST_domain_origin + rod2_r1 * ST_domain_len

    electrode_type = "rod"
    call CFG_add_get(cfg, "field_electrode_type", electrode_type, &
         "Electrode: sphere, rod, rod_cone_top, rod_rod, sphere_rod, user")
    !< [electrode_settings]

    if (associated(user_potential_bc)) then
       mg%sides_bc => user_potential_bc
    else
       ! Use one of the predefined boundary conditions
       select case (field_bc_type)
       case ("homogeneous")
          mg%sides_bc => field_bc_homogeneous
       case ("neumann")
          mg%sides_bc => field_bc_neumann
       case ("all_neumann")
          mg%sides_bc => field_bc_all_neumann
       case default
          error stop "field_bc_select error: invalid condition"
       end select
    end if

    ! Set the multigrid options. First define the variables to use
    mg%i_phi = i_phi
    mg%i_tmp = i_tmp
    mg%i_rhs = i_rhs

    if (ST_use_dielectric) tree%mg_i_eps = i_eps

    if (ST_use_electrode) then
       select case (electrode_type)
       case ("sphere")
          ! A single spherical electrode
          if (any(rod_r0 <= -1.0e10_dp)) &
               error stop "field_rod_r0 not set correctly"
          if (rod_radius <= 0) &
               error stop "field_rod_radius not set correctly"
          mg%lsf => sphere_lsf
       case ("rod")
          ! A single rod electrode with a semi-spherical cap
          call check_general_electrode_parameters()
          mg%lsf => rod_lsf
       case ("rod_cone_top")
          ! A single rod-shaped electrode with a conical top
          call check_general_electrode_parameters()
          if (cone_tip_radius <= 0 .or. cone_tip_radius > rod_radius) &
               error stop "cone_tip_radius should be smaller than rod radius"
          if (cone_length_frac < 0 .or. cone_length_frac > 1) &
               error stop "cone_length_frac not set correctly"

          call get_conical_rod_properties(rod_r0, rod_r1, rod_radius, &
               cone_tip_radius, cone_tip_center, cone_tip_r_curvature)

          mg%lsf => conical_rod_lsf
       case ("rod_rod")
          ! Two rod electrodes with semi-spherical caps
          call check_general_electrode_parameters()

          if (any(rod2_r0 <= -1.0e10_dp)) &
               error stop "field_rod2_r0 not set correctly"
          if (any(rod2_r1 <= -1.0e10_dp)) &
               error stop "field_rod2_r1 not set correctly"
          if (rod2_radius <= 0) &
               error stop "field_rod2_radius not set correctly"

          mg%lsf => rod_rod_lsf

          ! Provide a function to set the voltage on the electrodes
          mg%lsf_boundary_function => rod_rod_get_potential
       case ("sphere_rod")
          ! Sphere and rod electrode with semi-spherical cap
          call check_general_electrode_parameters()

          if (any(rod2_r0 <= -1.0e10_dp)) &
               error stop "field_rod2_r0 not set correctly"
          if (any(rod2_r1 <= -1.0e10_dp)) &
               error stop "field_rod2_r1 not set correctly"
          if (rod2_radius <= 0) &
               error stop "field_rod2_radius not set correctly"

          mg%lsf => sphere_rod_lsf
          mg%lsf_boundary_function => sphere_rod_get_potential
       case ("two_rod_cone_electrodes")
          ! Two rod-shaped electrodes with conical tops (for now assumed to have
          ! the same shape)
          call check_general_electrode_parameters()
          if (any(rod2_r0 <= -1.0e10_dp)) &
               error stop "field_rod2_r0 not set correctly"
          if (any(rod2_r1 <= -1.0e10_dp)) &
               error stop "field_rod2_r1 not set correctly"
          if (rod2_radius <= 0) &
               error stop "field_rod2_radius not set correctly"
          if (cone_tip_radius <= 0 .or. cone_tip_radius > rod_radius) &
               error stop "cone tip radius should be smaller than rod radius"
          if (cone2_tip_radius <= 0 .or. cone2_tip_radius > rod2_radius) &
               error stop "cone2 tip radius should be smaller than rod2 radius"
          if (cone_length_frac < 0 .or. cone_length_frac > 1) &
               error stop "cone_length_frac not set correctly"
          if (cone2_length_frac < 0 .or. cone2_length_frac > 1) &
               error stop "cone2_length_frac not set correctly"

          call get_conical_rod_properties(rod_r0, rod_r1, rod_radius, &
               cone_tip_radius, cone_tip_center, cone_tip_r_curvature)
          call get_conical_rod_properties(rod2_r0, rod2_r1, rod2_radius, &
               cone2_tip_radius, cone2_tip_center, cone2_tip_r_curvature)

          mg%lsf => two_conical_rods_lsf

          ! Provide a function to set the voltage on the electrodes
          mg%lsf_boundary_function => two_conical_rods_get_potential
       case ("coaxial")
          ! A coaxial geometry, with an inner and outer conductor

          if (rod_radius <= 0) &
               error stop "field_rod_radius not set correctly"
          if (rod2_radius <= 0) &
               error stop "field_rod2_radius not set correctly"

          ! Check if outer radius is correctly set
          if (ST_cylindrical) then
             if (rod2_radius > ST_domain_len(1)) &
                  error stop "field_rod2_radius too large"
          else if (NDIM == 2 .or. NDIM == 3) then
             if (rod2_radius > 0.5_dp * minval(ST_domain_len(1:2))) &
                  error stop "field_rod2_radius too large"
          end if

          if (NDIM == 3 .and. any(ST_periodic .neqv. [.false., .false., .true.])) &
               error stop "Coaxial requires periodic = F F T (last dimension periodic)"

          mg%lsf => coaxial_lsf
          mg%lsf_boundary_function => coaxial_get_potential

          ! Ensure that outer electrode is grounded
          field_electrode_grounded = .false.
          field_electrode2_grounded = .true.
          mg%sides_bc => field_bc_all_dirichlet
       case ("user")
          if (.not. associated(user_lsf)) then
             error stop "user_lsf not set"
          else
             mg%lsf => user_lsf
          end if

          if (associated(user_lsf_bc)) then
             mg%lsf_boundary_function => user_lsf_bc
          end if
       case default
          print *, "Electrode types: sphere, rod, rod_cone_top, rod_rod, user"
          error stop "Invalid electrode type"
       end select

       call af_set_cc_methods(tree, i_lsf, funcval=set_lsf_box)
       tree%mg_i_lsf = i_lsf

       mg%lsf_dist => mg_lsf_dist_gss

       if (rod_radius <= 0) then
          error stop "set field_rod_radius to smallest length scale of electrode"
       end if
       mg%lsf_length_scale = rod_radius
    end if

    call af_set_cc_methods(tree, i_electric_fld, &
         af_bc_neumann_zero, af_gc_interp)

  end subroutine field_initialize

  subroutine check_general_electrode_parameters()
    if (any(rod_r0 <= -1.0e10_dp)) &
         error stop "field_rod_r0 not set correctly"
    if (any(rod_r1 <= -1.0e10_dp)) &
         error stop "field_rod_r1 not set correctly"
    if (rod_radius <= 0) &
         error stop "field_rod_radius not set correctly"
  end subroutine check_general_electrode_parameters

  subroutine field_set_rhs(tree, s_in)
    use m_units_constants
    use m_chemistry
    use m_dielectric
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

    if (ST_use_dielectric) then
       call surface_surface_charge_to_rhs(tree, diel, i_surf_dens, i_rhs, fac)
    end if

  end subroutine field_set_rhs

  !> Compute electric field on the tree. First perform multigrid to get electric
  !> potential, then take numerical gradient to geld field.
  subroutine field_compute(tree, mg, s_in, time, have_guess)
    use m_chemistry
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg ! Multigrid option struct
    integer, intent(in)       :: s_in
    real(dp), intent(in)      :: time
    logical, intent(in)       :: have_guess
    integer                   :: i
    real(dp)                  :: max_rhs, residual_threshold, conv_fac
    real(dp)                  :: residual_ratio
    integer, parameter        :: max_initial_iterations = 100
    real(dp), parameter       :: max_residual = 1e8_dp
    real(dp), parameter       :: min_residual = 1e-6_dp
    real(dp)                  :: residuals(max_initial_iterations)

    call field_set_rhs(tree, s_in)
    call field_set_voltage(tree, time)

    call af_tree_maxabs_cc(tree, mg%i_rhs, max_rhs)

    ! With an electrode, the initial convergence testing should be less strict
    if (ST_use_electrode) then
       conv_fac = 1e-8_dp
    else
       conv_fac = 1e-10_dp
    end if

    ! Set threshold based on rhs and on estimate of round-off error, given by
    ! delta phi / dx^2 = (phi/L * dx)/dx^2
    ! Note that we use min_residual in case max_rhs and current_voltage are zero
    residual_threshold = max(min_residual, &
         max_rhs * ST_multigrid_max_rel_residual, &
         conv_fac * abs(current_voltage)/(ST_domain_len(NDIM) * af_min_dr(tree)))

    if (ST_use_electrode) then
       if (field_electrode_grounded) then
          mg%lsf_boundary_value = 0.0_dp
       else
          mg%lsf_boundary_value = current_voltage
       end if
    end if

    ! Perform a FMG cycle when we have no guess
    if (.not. have_guess) then
       do i = 1, max_initial_iterations
          call mg_fas_fmg(tree, mg, .true., .true.)
          call af_tree_maxabs_cc(tree, mg%i_tmp, residuals(i))

          if (residuals(i) < residual_threshold) then
             exit
          else if (i > 2) then
             ! Check if the residual is not changing much anymore, and if it is
             ! small enough, in which case we assume convergence
             residual_ratio = minval(residuals(i-2:i)) / &
                  maxval(residuals(i-2:i))
             if (residual_ratio < 2.0_dp .and. residual_ratio > 0.5_dp &
                  .and. residuals(i) < max_residual) exit
          end if
       end do

       ! Check for convergence
       if (i == max_initial_iterations + 1) then
          print *, "Iteration    residual"
          do i = 1, max_initial_iterations
             write(*, "(I4,E18.2)") i, residuals(i)
          end do
          print *, "Maybe increase multigrid_max_rel_residual?"
          error stop "No convergence in initial field computation"
       end if
    end if

    ! Perform cheaper V-cycles
    do i = 1, ST_multigrid_num_vcycles
       call mg_fas_vcycle(tree, mg, .true.)
       call af_tree_maxabs_cc(tree, mg%i_tmp, residuals(i))
       if (residuals(i) < residual_threshold) exit
    end do

    call field_from_potential(tree, mg)

  end subroutine field_compute

  !> Compute field from potential
  subroutine field_from_potential(tree, mg)
    use m_units_constants
    use m_dielectric
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg ! Multigrid option struct

    if (ST_use_dielectric) then
       call mg_compute_phi_gradient(tree, mg, electric_fld, -1.0_dp)
       call surface_correct_field_fc(tree, diel, i_surf_dens, &
            electric_fld, i_phi, UC_elem_charge / UC_eps0)
       call mg_compute_field_norm(tree, electric_fld, i_electric_fld)
    else
       call mg_compute_phi_gradient(tree, mg, electric_fld, -1.0_dp, i_electric_fld)
    end if

    ! Set the field norm also in ghost cells
    call af_gc_tree(tree, [i_electric_fld])
  end subroutine field_from_potential

  !> Set the current voltage
  subroutine field_set_voltage(tree, time)
    use m_units_constants
    use m_lookup_table
    use m_user_methods
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: electric_fld, t, tmp

    if (associated(user_field_amplitude)) then
       electric_fld = user_field_amplitude(tree, time)
       current_voltage = -ST_domain_len(NDIM) * electric_fld
       return
    end if

    select case (field_given_by)
    case (scalar_voltage)
       current_voltage = 0.0_dp

       if (time < field_pulse_period * field_num_pulses) then
          t = modulo(time, field_pulse_period)

          if (t < field_rise_time) then
             current_voltage = field_voltage * (t/field_rise_time)
          else if (t < field_pulse_width + field_rise_time) then
             current_voltage = field_voltage
          else
             tmp = t - (field_pulse_width + field_rise_time)
             current_voltage = field_voltage * max(0.0_dp, &
                  (1 - tmp/field_rise_time))
          end if
       end if
    case (tabulated_voltage)
       call LT_lin_interp_list(field_table_times, field_table_values, &
            time, current_voltage)
    end select
  end subroutine field_set_voltage

  !> Dirichlet boundary conditions for the potential in the last dimension,
  !> Neumann zero boundary conditions in the other directions
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
          bc_val  = current_voltage
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine field_bc_homogeneous

  !> A Dirichlet zero and non-zero Neumann boundary condition for the potential
  !> in the last dimension, Neumann zero boundary conditions in the other
  !> directions
  subroutine field_bc_neumann(box, nb, iv, coords, bc_val, bc_type)
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
          bc_type = af_bc_neumann
          bc_val  = current_voltage/ST_domain_len(NDIM)
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine field_bc_neumann

  !> All Neumann zero boundary conditions for the potential
  subroutine field_bc_all_neumann(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_neumann
    bc_val = 0.0_dp
  end subroutine field_bc_all_neumann

  !> All Dirichlet zero boundary conditions for the potential
  subroutine field_bc_all_dirichlet(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (ST_cylindrical) then
       if (nb == af_neighb_lowx) then
          bc_type = af_bc_neumann
       else
          bc_type = af_bc_dirichlet
       end if
       bc_val = 0.0_dp
    else
       bc_type = af_bc_dirichlet
       bc_val = 0.0_dp
    end if
  end subroutine field_bc_all_dirichlet

  ! This routine sets the level set function for a simple rod
  subroutine set_lsf_box(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = mg%lsf(rr)
    end do; CLOSE_DO
  end subroutine set_lsf_box

  real(dp) function sphere_lsf(r)
    real(dp), intent(in) :: r(NDIM)
    sphere_lsf = norm2(r - rod_r0) - rod_radius
  end function sphere_lsf

  real(dp) function rod_lsf(r)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)
    rod_lsf = GM_dist_line(r, rod_r0, rod_r1, NDIM) - rod_radius
  end function rod_lsf

  !> Compute several parameters for a conical rod
  subroutine get_conical_rod_properties(r0, r1, &
       rod_radius, tip_radius, cone_tip_center, cone_tip_r_curvature)
    real(dp), intent(in)  :: r0(NDIM) !< Beginning of rod
    real(dp), intent(in)  :: r1(NDIM) !< End of conical part
    real(dp), intent(in)  :: rod_radius   !< Radius of rod
    real(dp), intent(in)  :: tip_radius   !< Radius of curvature of tip
    real(dp), intent(out) :: cone_tip_center(NDIM)
    real(dp), intent(out) :: cone_tip_r_curvature
    real(dp)              :: cone_angle, cone_length

    ! Determine (half) the opening angle of the top cone, which goes from
    ! rod_radius to tip_radius over cone_length
    cone_length = cone_length_frac * norm2(r1 - r0)
    cone_angle  = atan((rod_radius - tip_radius) / cone_length)

    ! We have a point on a sphere with coordinates of the form (R*cos(a),
    ! R*sin(a)) = (tip_radius, y), so we can get R and subtract R sin(a) to
    ! obtain the center of the sphere
    cone_tip_r_curvature = tip_radius/cos(cone_angle)
    cone_tip_center = r1 - sin(cone_angle) * &
         cone_tip_r_curvature * (r1 - r0)/norm2(r1 - r0)
  end subroutine get_conical_rod_properties

  !> Helper function to compute lsf for a rod with a conical tip
  pure function conical_rod_lsf_arg(r, r0, r1, cone_tip_center, &
       rod_radius, tip_radius, cone_length_frac, r_curvature) result (lsf)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)
    real(dp), intent(in) :: r0(NDIM), r1(NDIM), cone_tip_center(NDIM)
    real(dp), intent(in) :: rod_radius, tip_radius, cone_length_frac
    real(dp), intent(in) :: r_curvature
    real(dp)             :: lsf
    real(dp)             :: dist_vec(NDIM), frac, radius_at_height, tmp

    ! Project onto line from r0 to r1
    call GM_dist_vec_line(r, r0, r1, NDIM, dist_vec, frac)

    if (frac <= 1 - cone_length_frac) then
       ! Cylindrical part
       lsf = norm2(dist_vec) - rod_radius
    else if (frac < 1.0_dp) then
       ! Conical part
       tmp = (1 - frac) / cone_length_frac ! between 0 and 1
       radius_at_height = tip_radius + tmp * (rod_radius - tip_radius)
       lsf = norm2(dist_vec) - radius_at_height
    else
       ! Spherical tip
       lsf = norm2(r - cone_tip_center) - r_curvature
    end if
  end function conical_rod_lsf_arg

  real(dp) function conical_rod_lsf(r)
    real(dp), intent(in)    :: r(NDIM)
    conical_rod_lsf = conical_rod_lsf_arg(r, rod_r0, rod_r1, cone_tip_center, &
         rod_radius, cone_tip_radius, cone_length_frac, cone_tip_r_curvature)
  end function conical_rod_lsf

  !> Get lsf for two conical rods
  real(dp) function two_conical_rods_lsf(r)
    real(dp), intent(in) :: r(NDIM)
    real(dp)             :: lsf_1, lsf_2

    lsf_1 = conical_rod_lsf_arg(r, rod_r0, rod_r1, cone_tip_center, &
         rod_radius, cone_tip_radius, cone_length_frac, cone_tip_r_curvature)
    lsf_2 = conical_rod_lsf_arg(r, rod2_r0, rod2_r1, cone2_tip_center, &
         rod2_radius, cone2_tip_radius, cone2_length_frac, cone2_tip_r_curvature)
    two_conical_rods_lsf = min(lsf_1, lsf_2)
  end function two_conical_rods_lsf

  real(dp) function two_conical_rods_get_potential(r) result(phi)
    real(dp), intent(in) :: r(NDIM)
    real(dp)             :: lsf_1, lsf_2

    lsf_1 = conical_rod_lsf_arg(r, rod_r0, rod_r1, cone_tip_center, &
         rod_radius, cone_tip_radius, cone_length_frac, cone_tip_r_curvature)
    lsf_2 = conical_rod_lsf_arg(r, rod2_r0, rod2_r1, cone2_tip_center, &
         rod2_radius, cone2_tip_radius, cone2_length_frac, cone2_tip_r_curvature)

    if (lsf_1 < lsf_2) then
       ! Closer to electrode 1
       if (field_electrode_grounded) then
          phi = 0.0_dp
       else
          phi = current_voltage
       end if
    else
       if (field_electrode2_grounded) then
          phi = 0.0_dp
       else
          phi = current_voltage
       end if
    end if
  end function two_conical_rods_get_potential

  !> Get level set function for case of two rods
  real(dp) function rod_rod_lsf(r)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)

    rod_rod_lsf = min(GM_dist_line(r, rod_r0, rod_r1, NDIM) - rod_radius, &
         GM_dist_line(r, rod2_r0, rod2_r1, NDIM) - rod2_radius)
  end function rod_rod_lsf

  !> Get potential to apply at electrode when there are two rods
  function rod_rod_get_potential(r) result(phi)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)
    real(dp)             :: phi, lsf_1, lsf_2

    ! Determine distance to electrodes
    lsf_1 = GM_dist_line(r, rod_r0, rod_r1, NDIM) - rod_radius
    lsf_2 = GM_dist_line(r, rod2_r0, rod2_r1, NDIM) - rod2_radius

    if (lsf_1 < lsf_2) then
       ! Closer to electrode 1
       if (field_electrode_grounded) then
          phi = 0.0_dp
       else
          phi = current_voltage
       end if
    else
       if (field_electrode2_grounded) then
          phi = 0.0_dp
       else
          phi = current_voltage
       end if
    end if
  end function rod_rod_get_potential

  !> Get level set function for case of sphere and rod
  real(dp) function sphere_rod_lsf(r)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)

    sphere_rod_lsf = min(sphere_lsf(r), &
         GM_dist_line(r, rod2_r0, rod2_r1, NDIM) - rod2_radius)
  end function sphere_rod_lsf

  !> Get potential to apply at electrode for sphere-rod case
  function sphere_rod_get_potential(r) result(phi)
    use m_geometry
    real(dp), intent(in) :: r(NDIM)
    real(dp)             :: phi, lsf_1, lsf_2

    ! Determine distance to electrodes
    lsf_1 = sphere_lsf(r)
    lsf_2 = GM_dist_line(r, rod2_r0, rod2_r1, NDIM) - rod2_radius

    if (lsf_1 < lsf_2) then
       ! Closer to electrode 1
       if (field_electrode_grounded) then
          phi = 0.0_dp
       else
          phi = current_voltage
       end if
    else
       if (field_electrode2_grounded) then
          phi = 0.0_dp
       else
          phi = current_voltage
       end if
    end if
  end function sphere_rod_get_potential

  subroutine coaxial_get_lsf(r, lsf_1, lsf_2)
    real(dp), intent(in)  :: r(NDIM)
    real(dp), intent(out) :: lsf_1, lsf_2
    real(dp)              :: domain_center(NDIM)

    if (ST_cylindrical) then
       domain_center(1) = 0
       domain_center(2) = ST_domain_origin(2) + 0.5_dp * ST_domain_len(2)
       lsf_1 = norm2(r(1:2) - domain_center(1:2)) - rod_radius
       lsf_2 = rod2_radius - norm2(r(1:2) - domain_center(1:2))
    else if (NDIM == 2 .or. NDIM == 3) then
       domain_center = ST_domain_origin + 0.5_dp * ST_domain_len
       ! Axis is parallel to last coordinate
       lsf_1 = norm2(r(1:2) - domain_center(1:2)) - rod_radius
       lsf_2 = rod2_radius - norm2(r(1:2) - domain_center(1:2))
    else
       error stop "coaxial not supported in 1D"
    end if
  end subroutine coaxial_get_lsf

  real(dp) function coaxial_lsf(r)
    real(dp), intent(in) :: r(NDIM)
    real(dp) :: lsf_1, lsf_2
    call coaxial_get_lsf(r, lsf_1, lsf_2)
    coaxial_lsf = min(lsf_1, lsf_2)
  end function coaxial_lsf

  !> Get potential to apply at electrode for sphere-rod case
  function coaxial_get_potential(r) result(phi)
    real(dp), intent(in) :: r(NDIM)
    real(dp)             :: phi, lsf_1, lsf_2

    call coaxial_get_lsf(r, lsf_1, lsf_2)

    if (lsf_1 < lsf_2) then
       phi = current_voltage
    else
       phi = 0.0_dp
    end if
  end function coaxial_get_potential

  !> Compute total field energy in Joule, defined as the volume integral over
  !> 1/2 * epsilon * E^2
  subroutine field_compute_energy(tree, field_energy)
    type(af_t), intent(in) :: tree
    real(dp), intent(out)  :: field_energy

    call af_reduction(tree, field_energy_box, reduce_sum, 0.0_dp, field_energy)
  end subroutine field_compute_energy

  !> Get the electrostatic field energy in a box
  real(dp) function field_energy_box(box)
    use m_units_constants
    type(box_t), intent(in) :: box
#if NDIM == 2
    integer                 :: i
    real(dp), parameter     :: twopi = 2 * acos(-1.0_dp)
#endif
    real(dp)                :: w(DTIMES(box%n_cell))
    integer                 :: nc

    nc = box%n_cell

    if (ST_use_dielectric) then
       w = 0.5_dp * UC_eps0 * box%cc(DTIMES(1:nc), i_eps) * product(box%dr)
    else
       w = 0.5_dp * UC_eps0 * product(box%dr)
    end if

#if NDIM == 2
    if (box%coord_t == af_cyl) then
       ! Weight by 2 * pi * r
       do i = 1, nc
          w(i, :) = w(i, :) * twopi * af_cyl_radius_cc(box, i)
       end do
    end if
#endif

    field_energy_box = sum(w * box%cc(DTIMES(1:nc), i_electric_fld)**2)
  end function field_energy_box

  real(dp) function reduce_sum(a, b)
    real(dp), intent(in) :: a, b
    reduce_sum = a + b
  end function reduce_sum

  function field_get_E_vector(box) result(E_vector)
    type(box_t), intent(in) :: box
    real(dp)                :: E_vector(DTIMES(1:box%n_cell), NDIM)
    integer                 :: nc

    nc = box%n_cell

#if NDIM == 1
    E_vector(DTIMES(:), 1) = 0.5_dp * (box%fc(1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1, electric_fld))
#elif NDIM == 2
    E_vector(DTIMES(:), 1) = 0.5_dp * (box%fc(1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1, electric_fld))
    E_vector(DTIMES(:), 2) = 0.5_dp * (box%fc(1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 2, electric_fld))
#elif NDIM == 3
    E_vector(DTIMES(:), 1) = 0.5_dp * (box%fc(1:nc, 1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1:nc, 1, electric_fld))
    E_vector(DTIMES(:), 2) = 0.5_dp * (box%fc(1:nc, 1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 1:nc, 2, electric_fld))
    E_vector(DTIMES(:), 3) = 0.5_dp * (box%fc(1:nc, 1:nc, 1:nc, 3, electric_fld) + &
         box%fc(1:nc, 1:nc, 2:nc+1, 3, electric_fld))
#endif
  end function field_get_E_vector

end module m_field
