#include "../afivo/src/cpp_macros.h"
!> Module to compute electric fields
module m_field
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
  real(dp) :: field_pulse_width = huge(1.0_dp)

  !> Number of voltage pulses
  integer :: field_num_pulses = 1

  !> Time of one complete voltage pulse
  real(dp) :: field_pulse_period = huge(1.0_dp)

  !> The (initial) vertical applied electric field
  real(dp) :: field_amplitude = undefined_real

  !> The applied voltage (vertical direction)
  real(dp) :: field_voltage = undefined_real

  !> Whether the electrode is grounded or at the applied voltage
  logical :: field_electrode_grounded = .false.

  !> First coordinate of electrode
  real(dp) :: field_rod_r0(NDIM) = -1.0_dp

  !> Second coordinate of electrode
  real(dp) :: field_rod_r1(NDIM) = -1.0_dp

  !> Third coordinate of electrode
  real(dp) :: field_rod_r2(NDIM) = -1.0_dp

  !> Electrode radius (in m, for standard rod electrode)
  real(dp) :: field_rod_radius = -1.0_dp

  !> Electrode tip radius (for conical electrode)
  real(dp) :: field_tip_radius = -1.0_dp

  ! Internal variables

  !> The current applied voltage
  real(dp), public, protected :: current_voltage = 0.0_dp

  ! Parameters that are pre-computed for conical rods
  ! @todo Add explanations/documentation for conical rods
  real(dp) :: conical_rod_R_o
  real(dp) :: conical_rod_y_o
  real(dp) :: conical_rod_ro(NDIM)

  character(string_len) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_set_rhs
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
    character(len=string_len)  :: given_by, user_value
    character(len=string_len)  :: electrode_type
    integer                    :: first_blank
    real(dp)                   :: R_rod, R_tip, y_tip, alpha
    real(dp)                   :: r0(NDIM), r1(NDIM), r2(NDIM)

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
       print *, "field_given_by value: ", trim(given_by)
       print *, "Options are: voltage, field, voltage_table, field_table"
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

    if (field_pulse_width < huge(1.0_dp) .and. field_rise_time <= 0) &
         error stop "Set field_rise_time when using field_pulse_width"

    if (field_num_pulses > 1) then
       if (field_pulse_period >= huge(1.0_dp)) &
            error stop "field_num_pulses > 1 requires field_pulse_period"
       if (field_pulse_width >= huge(1.0_dp)) &
            error stop "field_num_pulses > 1 requires field_pulse_width"
       if (field_pulse_width + 2 * field_rise_time > field_pulse_period) &
            error stop "field_pulse_period shorter than one full pulse"
    end if

    call CFG_add_get(cfg, "field_bc_type", field_bc_type, &
         "Boundary condition for electric potential")

    !< [electrode_settings]
    call CFG_add_get(cfg, "field_electrode_grounded", field_electrode_grounded, &
         "Whether the electrode is grounded or at the applied voltage")
    call CFG_add_get(cfg, "field_rod_r0", field_rod_r0, &
         "First electrode relative position")
    call CFG_add_get(cfg, "field_rod_r1", field_rod_r1, &
         "Second electrode relative position")
    call CFG_add_get(cfg, "field_rod_r2", field_rod_r2, &
         "Third electrode relative position")
    call CFG_add_get(cfg, "field_rod_radius", field_rod_radius, &
         "Electrode radius (in m, for standard rod electrode)")
    call CFG_add_get(cfg, "field_tip_radius", field_tip_radius, &
         "Electrode tip radius (for conical electrode)")

    electrode_type = "rod"
    call CFG_add_get(cfg, "field_electrode_type", electrode_type, &
         "Type of electrode (sphere, rod, rod_cone_top)")
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
          mg%lsf => sphere_lsf
       case ("rod")
          call check_rod_parameters()
          mg%lsf => rod_lsf
       case ("rod_cone_top")
          call check_rod_parameters()
          if (any(field_rod_r2 < field_rod_r1)) &
               error stop "Should have field_rod_r2 >= field_rod_r1"
          if (field_tip_radius <= 0) &
               error stop "field_tip_radius not set correctly"
          mg%lsf => conical_rod_top_lsf

          r0    = field_rod_r0 * ST_domain_len
          r1    = field_rod_r1 * ST_domain_len
          r2    = field_rod_r2 * ST_domain_len
          R_rod = field_rod_radius
          R_tip = field_tip_radius
          y_tip = (r0(NDIM)*R_rod-r1(NDIM)*R_tip)/(R_rod-R_tip)
          alpha = atan(R_rod/(r1(NDIM)-y_tip))

          conical_rod_R_o = R_tip/cos(alpha)
          conical_rod_y_o = r0(NDIM) + R_tip * tan(alpha)
          conical_rod_ro = [r0(1:NDIM-1), conical_rod_y_o]
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
          print *, "Electrode types: sphere, rod, rod_cone_top, user"
          error stop "Invalid electrode type"
       end select

       call af_set_cc_methods(tree, i_lsf, funcval=set_lsf_box)
       tree%mg_i_lsf = i_lsf

       mg%lsf_dist => mg_lsf_dist_gss

       if (field_rod_radius <= 0) then
          error stop "set field_rod_radius to smallest length scale of electrode"
       end if
       mg%lsf_length_scale = field_rod_radius
    end if

    call af_set_cc_methods(tree, i_electric_fld, &
         af_bc_neumann_zero, af_gc_interp)

  end subroutine field_initialize

  subroutine check_rod_parameters()
    if (any(field_rod_r0 <= -1.0_dp)) &
         error stop "field_rod_r0 not set correctly"
    if (any(field_rod_r1 <= -1.0_dp)) &
         error stop "field_rod_r1 not set correctly"
    if (field_rod_radius <= 0) &
         error stop "field_rod_radius not set correctly"
    if (any(field_rod_r1 < field_rod_r0)) &
         error stop "Should have field_rod_r1 >= field_rod_r0"
  end subroutine check_rod_parameters

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
    use m_units_constants
    use m_chemistry
    use m_dielectric
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

    ! Compute field from potential
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
  end subroutine field_compute

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
    sphere_lsf = norm2(r - field_rod_r0 * ST_domain_len) - field_rod_radius
  end function sphere_lsf

  real(dp) function rod_lsf(r)
    use m_geometry
    real(dp), intent(in)    :: r(NDIM)
    rod_lsf = GM_dist_line(r, field_rod_r0 * ST_domain_len, &
         field_rod_r1 * ST_domain_len, NDIM) - field_rod_radius
  end function rod_lsf

  function conical_rod_top_lsf(r) result(lsf)
    use m_geometry
    real(dp), intent(in)    :: r(NDIM)
    real(dp)                :: lsf
    real(dp)                :: r0(NDIM), r1(NDIM), r2(NDIM)
    real(dp)                :: radius_at_height

    r0 = field_rod_r0 * ST_domain_len
    r1 = field_rod_r1 * ST_domain_len
    r2 = field_rod_r2 * ST_domain_len

    if (r(NDIM) > r1(NDIM)) then
       lsf = GM_dist_line(r, r1, r2, NDIM) - field_rod_radius
    else if (r(NDIM) > r0(NDIM)) then
       radius_at_height = field_tip_radius + (r(NDIM) - r0(NDIM)) / &
            (r1(NDIM) - r0(NDIM)) * (field_rod_radius - field_tip_radius)
       lsf = GM_dist_line(r, r0, r1, NDIM) - radius_at_height
    else
       lsf = GM_dist_line(r, conical_rod_ro, conical_rod_ro, NDIM) - &
            conical_rod_R_o
    end if
  end function conical_rod_top_lsf

end module m_field
