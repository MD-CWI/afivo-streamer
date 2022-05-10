#include "../afivo/src/cpp_macros.h"
!> Module to compute electric fields
module m_field
  use m_af_all
  use m_types
  use m_streamer

  implicit none
  private


  !> Use a table with time vs voltage
  logical :: voltage_table_use

  !> List of times
  real(dp), allocatable :: voltage_table_times(:)

  !> List of voltages
  real(dp), allocatable :: voltage_table_voltages(:)

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

  !> Amplitude of sinusoidal modification
  real(dp) :: field_sin_amplitude = 0.0_dp

  !> Frequency (Hz) of sinusoidal modification
  real(dp) :: field_sin_freq = 0.0_dp

  !> Linear derivative of background field
  real(dp) :: field_lin_deriv = 0.0_dp

  !> Decay time of background field
  real(dp) :: field_decay_time = huge(1.0_dp)

  !> Linear rise time of field (s)
  real(dp) :: field_rise_time = 0.0_dp

  !> The (initial) vertical applied electric field
  real(dp), public, protected :: field_amplitude = 1.0e6_dp

  !> The current vertical applied electric field
  real(dp), public, protected :: current_field_amplitude

  !> The applied voltage (vertical direction)
  real(dp), public, protected :: field_voltage

  !> Whether the voltage has been set externally
  logical :: voltage_set_externally = .false.

  !> Whether the electrode is grounded or at the applied voltage
  logical, public, protected :: field_electrode_grounded = .false.

  !> First electrode position
  real(dp), public, protected :: field_rod_r0(NDIM) = -1.0_dp

  !> Second electrode position
  real(dp), public, protected :: field_rod_r1(NDIM) = -1.0_dp

  !> Third electrode position
  real(dp), public, protected :: field_rod_r2(NDIM) = -1.0_dp

  !> Electrode radius (in m, for standard rod electrode)
  real(dp), public, protected :: field_rod_radius = -1.0_dp

  !> Electrode voltage (in V)
  real(dp), public, protected :: electrode_voltage = 0.0_dp

  !> Whether the voltage of the electrode is used or the current field amplitude
  logical, public, protected :: use_electrode_voltage = .false.

  !> Rise/fall time of electrode voltage
  real(dp), public, protected :: electrode_voltage_t_risefall = -1

  !> Time of constant applied voltage within a pulse
  real(dp), public, protected :: electrode_voltage_t_pulse = -1
  
  !> Time between 2 voltage pulses
  real(dp), public, protected :: electrode_voltage_t_inter = -1

  !> Electrode tip radius (for conical electrode)
  real(dp), public, protected :: field_tip_radius = -1.0_dp

  logical  :: field_stability_search    = .false.
  real(dp) :: field_stability_zmin      = 0.2_dp
  real(dp) :: field_stability_zmax      = 1.0_dp
  real(dp) :: field_stability_threshold = 3e6_dp

  real(dp) :: field_point_charge = 0.0_dp
  real(dp) :: field_point_r0(NDIM) = 0.0_dp

  ! Parameters that are pre-computed for conical rods
  ! @todo Add explanations/documentation for conical rods
  real(dp) :: conical_rod_R_o
  real(dp) :: conical_rod_y_o
  real(dp) :: conical_rod_ro(NDIM)

  character(string_len) :: field_bc_type = "homogeneous"

  public :: field_initialize
  public :: field_compute
  public :: field_set_rhs
  public :: field_get_amplitude
  public :: field_set_voltage
  public :: field_set_voltage_externally

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
    character(len=string_len)  :: field_table, electrode_type, voltage_table
    real(dp)                   :: R_rod, R_tip, y_tip, alpha
    real(dp)                   :: r0(NDIM), r1(NDIM), r2(NDIM)

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

    voltage_table = undefined_str
    call CFG_add_get(cfg, "voltage_table", voltage_table, &
         "File containing time vs voltage (V)")
    if (voltage_table /= undefined_str) then
      voltage_table_use = .true.
      call table_from_file(voltage_table, "time_vs_voltage", &
            voltage_table_times, voltage_table_voltages)
    else
      voltage_table_use = .false.
    end if

    call CFG_add_get(cfg, "field_mod_t0", field_mod_t0, &
         "Modify electric field after this time (s)")
    call CFG_add_get(cfg, "field_mod_t1", field_mod_t1, &
         "Modify electric field up to this time (s)")
    call CFG_add_get(cfg, "field_sin_amplitude", field_sin_amplitude, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add_get(cfg, "field_sin_freq", field_sin_freq, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add_get(cfg, "field_lin_deriv", field_lin_deriv, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add_get(cfg, "field_decay_time", field_decay_time, &
         "Decay time of field (s)")
    call CFG_add_get(cfg, "field_rise_time", field_rise_time, &
         "Linear rise time of field (s)")
    call CFG_add_get(cfg, "field_amplitude", field_amplitude, &
         "The (initial) vertical applied electric field (V/m)")
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
    field_point_r0(:) = 0.0_dp
    field_point_r0(NDIM) = -1.0_dp
    call CFG_add_get(cfg, "field_point_r0", field_point_r0, &
         "Relative position of point charge (outside domain)")

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
    call CFG_add_get(cfg, "electrode_voltage", electrode_voltage, &
         "The (initial) applied electrode voltage (V).")
    call CFG_add_get(cfg, "use_electrode_voltage", use_electrode_voltage, &
         "Whether the defined electrode voltage is used instead of the field amplitude.")
    call CFG_add_get(cfg, "electrode_voltage_t_risefall", electrode_voltage_t_risefall, &
         "Rise/fall time of electrode voltage")
    call CFG_add_get(cfg, "electrode_voltage_t_pulse", electrode_voltage_t_pulse, &
         "Time of constant applied voltage within a pulse")
    call CFG_add_get(cfg, "electrode_voltage_t_inter", electrode_voltage_t_inter, &
         "Time between 2 voltage pulses")
    call CFG_add_get(cfg, "field_tip_radius", field_tip_radius, &
         "Electrode tip radius (for conical electrode)")

    electrode_type = "rod"
    call CFG_add_get(cfg, "field_electrode_type", electrode_type, &
         "Type of electrode")

    if (associated(user_potential_bc)) then
       mg%sides_bc => user_potential_bc
    else
       ! Use one of the predefined boundary conditions
       select case (field_bc_type)
       case ("homogeneous")
          mg%sides_bc => field_bc_homogeneous
       case ("neumann")
          mg%sides_bc => field_bc_neumann
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

    if (ST_use_dielectric) tree%mg_i_eps = i_eps

    if (ST_use_electrode) then
       if (any(field_rod_r0 <= -1.0_dp)) &
            error stop "field_rod_r0 not set correctly"
       if (any(field_rod_r1 <= -1.0_dp)) &
            error stop "field_rod_r1 not set correctly"
       if (field_rod_radius <= 0) &
            error stop "field_rod_radius not set correctly"
       if (any(field_rod_r1 < field_rod_r0)) &
            error stop "Should have field_rod_r1 >= field_rod_r0"

       select case (electrode_type)
       case ("rod")
          call af_set_cc_methods(tree, i_lsf, funcval=set_field_rod_lsf)
          mg%lsf => field_rod_lsf
       case ("rod_cone_top")
          if (any(field_rod_r2 < field_rod_r1)) &
               error stop "Should have field_rod_r2 >= field_rod_r1"
          if (field_tip_radius <= 0) &
               error stop "field_tip_radius not set correctly"
          call af_set_cc_methods(tree, i_lsf, funcval=set_field_conical_rod_top_lsf)
          mg%lsf => field_conical_rod_top_lsf

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
       case default
          error stop "Invalid electrode type (option: rod, rod_cone_top)"
       end select

       tree%mg_i_lsf = i_lsf
    end if

    mg%lsf_dist => mg_lsf_dist_gss
    mg%lsf_length_scale = field_rod_radius

    call af_set_cc_methods(tree, i_electric_fld, &
         af_bc_neumann_zero, af_gc_interp)

  end subroutine field_initialize

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
    ! Note that we use min_residual in case max_rhs and field_voltage are zero
    residual_threshold = max(min_residual, &
         max_rhs * ST_multigrid_max_rel_residual, &
         conv_fac * abs(field_voltage)/(ST_domain_len(NDIM) * af_min_dr(tree)))

    if (ST_use_electrode) then
       if (field_electrode_grounded) then
          mg%lsf_boundary_value = 0.0_dp
       else
          mg%lsf_boundary_value = field_voltage
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

    select case (ST_field_correction)
    case ("divE")
       call af_loop_box(tree, correct_field_divE_box)
    case ("harmonic")
       call af_loop_box(tree, correct_field_harmonic_box)
    case ("none")
       continue
    case default
       error stop "Unknown fixes%field_correction"
    end select

    ! Set the field norm also in ghost cells
    call af_gc_tree(tree, [i_electric_fld])
  end subroutine field_compute

  subroutine correct_field_divE_box(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK
    real(dp)                   :: tmp, Elo(NDIM), Ehi(NDIM)
    real(dp)                   :: Emin(NDIM), dE(NDIM)
    real(dp)                   :: eps = 1e-10_dp
    integer                    :: nc

    nc = box%n_cell

    do KJI_DO(1, nc)
       ! Compute f = 1 - |div(E)| / (|d/dx Ex| + |d/dy Ey|)
       ! Then take min(E) + f * delta_E (component wise)
       ! Todo: think about unequal mesh spacing

#if NDIM == 2
       ! Fields on lower and upper cell faces
       Elo = box%fc(IJK, 1:NDIM, electric_fld)
       Ehi = [box%fc(i+1, j, 1, electric_fld), &
            box%fc(i, j+1, 2, electric_fld)]
#endif

       ! Emin contains the field components of smallest amplitude
       where (abs(Elo) < abs(Ehi))
          Emin = Elo
          dE = Ehi - Elo
       elsewhere
          Emin = Ehi
          dE = Elo - Ehi
       end where

       tmp = 1 - abs(sum(Ehi-Elo)) / (eps + sum(abs(dE)))
       box%cc(IJK, i_electric_fld) = norm2(Emin + 0.5_dp * tmp * dE)
    end do; CLOSE_DO
  end subroutine correct_field_divE_box

  subroutine correct_field_harmonic_box(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK
    real(dp)                   :: f(2*NDIM), hmean
    real(dp)                   :: eps = 1e-10_dp
    integer                    :: nc

    nc = box%n_cell

    do KJI_DO(1, nc)
#if NDIM == 1
       f = abs([box%fc(i, 1, electric_fld), &
            box%fc(i+1, 1, electric_fld)])
       hmean = 2*f(1)*f(2)/(f(1)+f(2)+eps)
#elif NDIM == 2
       f = abs([box%fc(i, j, 1, electric_fld), &
            box%fc(i+1, j, 1, electric_fld), &
            box%fc(i, j, 2, electric_fld), &
            box%fc(i, j+1, 2, electric_fld)])
       hmean = norm2([2*f(1)*f(2)/(f(1)+f(2)+eps), &
            2*f(3)*f(4)/(f(3)+f(4)+eps)])
#elif NDIM == 3
       f = abs([&
            box%fc(i, j, k, 1, electric_fld), &
            box%fc(i+1, j, k, 1, electric_fld), &
            box%fc(i, j, k, 2, electric_fld), &
            box%fc(i, j+1, k, 2, electric_fld), &
            box%fc(i, j, k, 3, electric_fld), &
            box%fc(i, j, k+1, 3, electric_fld)])
       hmean = norm2([&
            2*f(1)*f(2)/(f(1)+f(2)+eps), &
            2*f(3)*f(4)/(f(3)+f(4)+eps), &
            2*f(5)*f(6)/(f(5)+f(6)+eps)])
#endif
       box%cc(IJK, i_electric_fld) = hmean
    end do; CLOSE_DO
  end subroutine correct_field_harmonic_box

  !> Compute the electric field at a given time
  function field_get_amplitude(tree, time) result(electric_fld)
    use m_units_constants
    use m_lookup_table
    use m_user_methods
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: electric_fld, t_rel

    if (associated(user_field_amplitude)) then
       electric_fld = user_field_amplitude(tree, time)
    else if (field_table_use) then
       call LT_lin_interp_list(field_table_times, field_table_fields, &
            time, electric_fld)
    else
       electric_fld = field_amplitude

       if (time < field_rise_time) then
          electric_fld = electric_fld * (time/field_rise_time)
       end if

       ! TODO: simplify stuff below
       t_rel = time - field_mod_t0
       t_rel = min(t_rel, field_mod_t1-field_mod_t0)

       if (t_rel > 0) then
          electric_fld = electric_fld * exp(-t_rel/field_decay_time) + &
               t_rel * field_lin_deriv + &
               field_sin_amplitude * &
               sin(t_rel * field_sin_freq * 2 * UC_pi)
       end if
    end if

  end function field_get_amplitude

  !> Compute the voltage at a given time
  subroutine field_set_voltage(tree, time)
   use m_lookup_table

    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: electrode_voltage_t_block
    real(dp)               :: electrode_voltage_period
    real(dp)               :: time_voltage_shifted

    if (.not. voltage_set_externally) then
       current_field_amplitude = field_get_amplitude(tree, time)
       if (.not. ST_use_electrode) then
         field_voltage = -ST_domain_len(NDIM) * current_field_amplitude
       else

         if (voltage_table_use) then
            call LT_lin_interp_list(voltage_table_times, voltage_table_voltages, &
               time, field_voltage)
         else if (use_electrode_voltage) then

            if (electrode_voltage_t_inter >= 0 .and. electrode_voltage_t_pulse >= 0 .and. electrode_voltage_t_risefall >= 0) then
               !> We will continuously generate block pulses with a rise and fall time.
               electrode_voltage_t_block = 2 * electrode_voltage_t_risefall + electrode_voltage_t_pulse
               electrode_voltage_period = electrode_voltage_t_block + electrode_voltage_t_inter
               
               !> We calculate where our current time is in the range [0, electrode_voltage_period]
               time_voltage_shifted = modulo(time, electrode_voltage_period)

                !> Rise of voltage
               if (time_voltage_shifted <= electrode_voltage_t_risefall) then
                  field_voltage = electrode_voltage * (time_voltage_shifted / electrode_voltage_t_risefall)
               !> Constant voltage within a voltage block
               else if (time_voltage_shifted <= (electrode_voltage_t_risefall + electrode_voltage_t_pulse)) then
                  field_voltage = electrode_voltage
               !> Fall of voltage
               else if (time_voltage_shifted <= electrode_voltage_t_block) then
                  field_voltage = electrode_voltage * &
                     (1 - (time_voltage_shifted - (electrode_voltage_t_pulse + electrode_voltage_t_risefall)) &
                                                                              / electrode_voltage_t_risefall)
               !> Inter voltage block 0 V time
               else if (time_voltage_shifted <= (electrode_voltage_t_block + electrode_voltage_t_inter)) then
                  field_voltage = 0
               end if

            else
               field_voltage = electrode_voltage
            end if
         else
            field_voltage = -(ST_domain_len(NDIM) * (1 - (abs(field_rod_r0(NDIM) - field_rod_r1(NDIM)))))&
                         * current_field_amplitude
         end if
      end if
    end if
  end subroutine field_set_voltage

  !> Set the voltage
  subroutine field_set_voltage_externally(voltage)
    real(dp), intent(in) :: voltage
    voltage_set_externally = .true.
    field_voltage = voltage

   if (.not. ST_use_electrode) then
      current_field_amplitude = -voltage / ST_domain_len(NDIM)
   else
      current_field_amplitude = -voltage / (ST_domain_len(NDIM) * (1 - (abs(field_rod_r0(NDIM) - field_rod_r1(NDIM)))))
   end if
  end subroutine field_set_voltage_externally

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
          bc_val  = field_voltage
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
          bc_val  = -current_field_amplitude
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine field_bc_neumann

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

  ! This routine sets the level set function for a simple rod
  subroutine set_field_rod_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = field_rod_lsf(rr)
    end do; CLOSE_DO

  end subroutine set_field_rod_lsf

  real(dp) function field_rod_lsf(r)
    use m_geometry
    real(dp), intent(in)    :: r(NDIM)
    field_rod_lsf = GM_dist_line(r, field_rod_r0 * ST_domain_len, &
         field_rod_r1 * ST_domain_len, NDIM) - field_rod_radius
  end function field_rod_lsf

  ! This routine sets the level set function for a simple rod
  subroutine set_field_conical_rod_top_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = field_conical_rod_top_lsf(rr)
    end do; CLOSE_DO

  end subroutine set_field_conical_rod_top_lsf

  function field_conical_rod_top_lsf(r) result(lsf)
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
  end function field_conical_rod_top_lsf

end module m_field
