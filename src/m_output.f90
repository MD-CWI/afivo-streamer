#include "../afivo/src/cpp_macros.h"
!> Module with methods and parameters related to writing output
module m_output
  use m_af_all
  use m_types
  use m_streamer
  use m_gas
  use m_lookup_table
  use m_units_constants

  implicit none
  private

  ! Name for the output files
  character(len=string_len), public, protected :: output_name = "output/my_sim"

  ! Optional variable (to show ionization source term)
  integer, public, protected :: i_src = -1 ! Source term

  ! If defined, only output these variables
  character(len=af_nlen), allocatable :: output_only(:)

  ! Time between writing output
  real(dp), public, protected :: output_dt = 1.0e-10_dp

  ! Include ionization source term in output
  logical, public, protected :: output_src_term = .false.

  ! If positive: decay rate for source term (1/s) for time-averaged values
  real(dp), public, protected :: output_src_decay_rate = -1.0_dp

  ! Write to a log file for regression testing
  logical, protected :: output_regression_test = .false.

  ! Write output along a line
  logical, public, protected :: lineout_write = .false.

  ! Write Silo output
  logical, public, protected :: silo_write = .true.

  ! Write Silo output every N outputs
  integer, public, protected :: silo_per_outputs = 1

  ! Write binary output files (to resume later)
  logical, public, protected :: datfile_write = .false.

  ! Write binary output files every N outputs
  integer, public, protected :: datfile_per_outputs = 1

  ! Use this many points for lineout data
  integer, public, protected :: lineout_npoints = 500

  ! Relative position of line minimum coordinate
  real(dp), public, protected :: lineout_rmin(3) = 0.0_dp

  ! Relative position of line maximum coordinate
  real(dp), public, protected :: lineout_rmax(3) = 1.0_dp

  ! Which variable to include in lineout
  integer, allocatable, public, protected :: lineout_ivar(:)

  ! Write integral over cross-section data output
  logical, public, protected :: cross_write = .false.

  ! Integrate up to this r value
  real(dp), public, protected :: cross_rmax = 2.0e-3_dp

  ! Use this many points for cross-section data
  integer, public, protected :: cross_npoints = 500

  ! Write uniform output in a plane
  logical, public, protected :: plane_write = .false.

  ! Which variable to include in plane
  integer, allocatable, public, protected :: plane_ivar(:)

  ! Use this many points for plane data
  integer, public, protected :: plane_npixels(2) = [64, 64]

  ! Relative position of plane minimum coordinate
  real(dp), public, protected :: plane_rmin(NDIM) = 0.0_dp

  ! Relative position of plane maximum coordinate
  real(dp), public, protected :: plane_rmax(NDIM) = 1.0_dp

  ! Output electric field maxima and their locations
  logical, public, protected :: field_maxima_write = .false.

  ! Threshold value (V/m) for electric field maxima
  real(dp), public, protected :: field_maxima_threshold = 0.0_dp

  ! Minimal distance (m) between electric field maxima
  real(dp), public, protected :: field_maxima_distance = 0.0_dp

  ! Print status every this many seconds
  real(dp), public, protected :: output_status_delay = 60.0_dp

  ! To reduce output when the voltage is off
  integer, public, protected :: output_dt_factor_pulse_off = 1

  ! Density threshold for detecting plasma regions
  real(dp) :: density_threshold = 1e18_dp

  ! Number of extra variables to add to output
  integer :: n_extra_vars = 0

  ! Names of extra variables to add to output
  character(len=af_nlen) :: output_extra_vars(100)

  ! Output the electron energy in eV from the local field approximation
  logical :: output_electron_energy = .false.

  ! Output the conductivity of the plasma
  logical :: output_conductivity = .false.

  ! Table with energy in eV vs electric field
  type(LT_t) :: eV_vs_fld

  ! Public methods
  public :: output_initialize
  public :: output_initial_summary
  public :: output_write
  public :: output_log
  public :: output_status

contains

  subroutine output_initialize(tree, cfg)
    use m_config
    use m_table_data
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    character(len=name_len), allocatable :: varname(:)
    character(len=af_nlen)     :: empty_names(0)
    integer                    :: n, i
    character(len=string_len)  :: td_file
    real(dp)                   :: tmp
    real(dp), allocatable      :: x_data(:), y_data(:)

    call CFG_add_get(cfg, "output%name", output_name, &
         "Name for the output files (e.g. output/my_sim)")

    call check_path_writable(output_name)

    call CFG_add_get(cfg, "output%status_delay", output_status_delay, &
         "Print status every this many seconds")

    call CFG_add_get(cfg, "output%dt_factor_pulse_off", &
         output_dt_factor_pulse_off, &
         "To reduce output when the voltage is off")

    call CFG_add(cfg, "output%only", empty_names, &
         "If defined, only output these variables", .true.)
    call CFG_get_size(cfg, "output%only", n)
    allocate(output_only(n))
    call CFG_get(cfg, "output%only", output_only)

    if (size(output_only) > 0) then
       tree%cc_write_output(:) = .false.
       do n = 1, size(output_only)
          i = af_find_cc_variable(tree, trim(output_only(n)))
          tree%cc_write_output(i) = .true.
       end do
    end if

    call CFG_add_get(cfg, "output%dt", output_dt, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "output%src_term", output_src_term, &
         "Include ionization source term in output")
    call CFG_add_get(cfg, "output%regression_test", output_regression_test, &
         "Write to a log file for regression testing")
    call CFG_add_get(cfg, "output%density_threshold", density_threshold, &
         "Electron density threshold (1/m3, will be scaled by N)")

    if (output_src_term) then
       call af_add_cc_variable(tree, "src", ix=i_src)
       call af_set_cc_methods(tree, i_src, af_bc_neumann_zero, af_gc_interp)
    end if

    tmp = 1/output_src_decay_rate
    call CFG_add_get(cfg, "output%src_decay_time", tmp, &
         "If positive: decay time for source term (s) for time-averaged values")
    output_src_decay_rate = 1/tmp

    call CFG_add_get(cfg, "silo_write", silo_write, &
         "Write silo output")
    call CFG_add_get(cfg, "silo%per_outputs", silo_per_outputs, &
         "Write silo output files every N outputs")

    call CFG_add_get(cfg, "datfile%write", datfile_write, &
         "Write binary output files (to resume later)")
    call CFG_add_get(cfg, "datfile%per_outputs", datfile_per_outputs, &
         "Write binary output files every N outputs")

    call CFG_add_get(cfg, "lineout%write", lineout_write, &
         "Write output along a line")
    if (lineout_write) then
       call CFG_add(cfg, "lineout%varname", ["e"], &
            "Names of variable to write in lineout", dynamic_size=.true.)
       call CFG_get_size(cfg, "lineout%varname", n)
       allocate(varname(n))
       call CFG_get(cfg, "lineout%varname", varname)

       allocate(lineout_ivar(n))
       do i = 1, n
          lineout_ivar(i) = af_find_cc_variable(tree, trim(varname(i)))
       end do
       deallocate(varname)

       call CFG_add_get(cfg, "lineout%npoints", lineout_npoints, &
            "Use this many points for lineout data")
       call CFG_add_get(cfg, "lineout%rmin", lineout_rmin(1:NDIM), &
            "Relative position of line minimum coordinate")
       call CFG_add_get(cfg, "lineout%rmax", lineout_rmax(1:NDIM), &
            "Relative position of line maximum coordinate")
    end if

    call CFG_add_get(cfg, "cross%write", cross_write, &
         "Write integral over cross-section data output")
    call CFG_add_get(cfg, "cross%rmax", cross_rmax, &
         "Integrate up to this r value")
    call CFG_add_get(cfg, "cross%npoints", cross_npoints, &
         "Use this many points for cross-section data")

    call CFG_add_get(cfg, "plane%write", plane_write, &
         "Write uniform output in a plane")

    if (plane_write) then
       call CFG_add(cfg, "plane%varname", ["e"], &
            "Names of variable to write in a plane", dynamic_size=.true.)

       call CFG_get_size(cfg, "plane%varname", n)
       allocate(varname(n))
       call CFG_get(cfg, "plane%varname", varname)

       allocate(plane_ivar(n))
       do i = 1, n
          plane_ivar(i) = af_find_cc_variable(tree, trim(varname(i)))
       end do
       deallocate(varname)

       call CFG_add_get(cfg, "plane%npixels", plane_npixels, &
            "Use this many pixels for plane data")
       call CFG_add_get(cfg, "plane%rmin", plane_rmin(1:NDIM), &
            "Relative position of plane minimum coordinate")
       call CFG_add_get(cfg, "plane%rmax", plane_rmax(1:NDIM), &
            "Relative position of plane maximum coordinate")
    end if

    call CFG_add_get(cfg, "field_maxima%write", field_maxima_write, &
         "Output electric field maxima and their locations")
    call CFG_add_get(cfg, "field_maxima%threshold", field_maxima_threshold, &
         "Threshold value (V/m) for electric field maxima")
    call CFG_add_get(cfg, "field_maxima%distance", field_maxima_distance, &
         "Minimal distance (m) between electric field maxima")

    call CFG_add_get(cfg, "output%electron_energy", output_electron_energy, &
         "Show the electron energy in eV from the local field approximation")
    call CFG_add_get(cfg, "output%conductivity", output_conductivity, &
         "Output the conductivity of the plasma")

    if (output_electron_energy) then
       n_extra_vars = n_extra_vars + 1
       output_extra_vars(n_extra_vars) = "eV"

       ! Create a lookup table for the model coefficients
       eV_vs_fld = LT_create(table_min_townsend, table_max_townsend, &
            table_size, 1, table_xspacing)
       call CFG_get(cfg, "input_data%file", td_file)

       ! Read table with E/N vs electron energy (eV)
       call table_from_file(td_file, "Mean energy (eV)", x_data, y_data)
       call table_set_column(eV_vs_fld, 1, x_data, y_data)
    end if

    if (output_conductivity) then
       n_extra_vars = n_extra_vars + 1
       output_extra_vars(n_extra_vars) = "sigma"
    end if

    call CFG_add(cfg, "output%write_derivatives", empty_names, &
         "Write derivatives of these species to output", dynamic_size=.true.)
    call CFG_get_size(cfg, "output%write_derivatives", n)
    allocate(varname(n))
    call CFG_get(cfg, "output%write_derivatives", varname)

    do i = 1, n
       n_extra_vars = n_extra_vars + 1
       output_extra_vars(n_extra_vars) = "src_" // trim(varname(i))
    end do
    deallocate(varname)

  end subroutine output_initialize

  !> Write a summary of the model and parameters used
  subroutine output_initial_summary(tree)
    use m_chemistry
    type(af_t), intent(in) :: tree
    character(len=string_len) :: fname

    fname = trim(output_name) // "_summary.txt"
    call chemistry_write_summary(trim(fname))
    call output_stoichiometric_matrix(trim(output_name)//"_stoich_matrix.txt")
    call output_chemical_species(trim(output_name)//"_species.txt")
    call output_chemical_reactions(trim(output_name)//"_reactions.txt")
    call output_chemical_rates(trim(output_name)//"_rates.txt", .true.)
    call output_chemical_amounts(tree, trim(output_name)//"_amounts.txt", .true.)
  end subroutine output_initial_summary

  subroutine output_write(tree, output_cnt, wc_time, write_sim_data)
    use m_photoi
    use m_field
    use m_user_methods
    use m_analysis
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: output_cnt
    real(dp), intent(in)      :: wc_time
    interface
       subroutine write_sim_data(my_unit)
         integer, intent(in) :: my_unit
       end subroutine write_sim_data
    end interface
    integer                   :: i
    character(len=string_len) :: fname

    if (compute_power_density) then
       call af_loop_box(tree, set_power_density_box)
    end if

    if (gas_dynamics) then
       call af_loop_box(tree, set_gas_primitives_box)
    end if

    if (silo_write .and. &
         modulo(output_cnt, silo_per_outputs) == 0) then
       ! Because the mesh could have changed
       if (photoi_enabled) call photoi_set_src(tree, global_dt)
       call field_set_rhs(tree, 0)

       do i = 1, tree%n_var_cell
          if (tree%cc_write_output(i) .and. .not. &
               associated(tree%cc_methods(i)%funcval)) then
             call af_restrict_tree(tree, [i])
             call af_gc_tree(tree, [i])
          end if
       end do

       write(fname, "(A,I6.6)") trim(output_name) // "_", output_cnt
       if (n_extra_vars > 0) then
          call af_write_silo(tree, fname, output_cnt, global_time, &
               add_vars=add_variables, &
               add_names=output_extra_vars(1:n_extra_vars))
       else
          call af_write_silo(tree, fname, output_cnt, global_time)
       end if
    end if

    if (datfile_write .and. &
         modulo(output_cnt, datfile_per_outputs) == 0) then
       call af_write_tree(tree, fname, write_sim_data)
    end if

    write(fname, "(A,I6.6)") trim(output_name) // "_log.txt"
    if (associated(user_write_log)) then
       call user_write_log(tree, fname, output_cnt)
    else
       call output_log(tree, fname, output_cnt, wc_time)
    end if

    call output_chemical_rates(trim(output_name)//"_rates.txt", .false.)
    call output_chemical_amounts(tree, trim(output_name)//"_amounts.txt", .false.)

    if (output_regression_test) then
       write(fname, "(A,I6.6)") trim(output_name) // "_rtest.log"
       call output_regression_log(tree, fname, output_cnt, wc_time)
    end if

    if (field_maxima_write) then
       write(fname, "(A,I6.6,A)") trim(output_name) // &
            "_Emax_", output_cnt, ".txt"
       call output_fld_maxima(tree, fname)
    end if

#if NDIM > 1
    if (plane_write) then
       write(fname, "(A,I6.6)") trim(output_name) // &
            "_plane_", output_cnt
       call af_write_plane(tree, fname, plane_ivar, &
            plane_rmin * ST_domain_len + ST_domain_origin, &
            plane_rmax * ST_domain_len + ST_domain_origin, &
            plane_npixels)
    end if
#endif

    if (lineout_write) then
       write(fname, "(A,I6.6)") trim(output_name) // &
            "_line_", output_cnt
       call af_write_line(tree, trim(fname), lineout_ivar, &
            r_min = lineout_rmin(1:NDIM) * ST_domain_len + ST_domain_origin, &
            r_max = lineout_rmax(1:NDIM) * ST_domain_len + ST_domain_origin, &
            n_points=lineout_npoints)
    end if

#if NDIM == 2
    if (ST_cylindrical) then
       if(cross_write) then
          write(fname, "(A,I6.6)") trim(output_name) // &
               "_cross_", output_cnt
          call output_cross(tree, fname)
       end if
    end if
#endif

  end subroutine output_write

  subroutine add_variables(box, new_vars, n_var)
    use m_gas
    use m_transport_data
    use m_chemistry
    type(box_t), intent(in) :: box
    integer, intent(in)     :: n_var
    real(dp)                :: new_vars(DTIMES(0:box%n_cell+1), n_var)
    integer                 :: n, nc, n_cells, i_species
    character(len=name_len) :: species_name
    real(dp)                :: N_inv(DTIMES(0:box%n_cell+1))
    real(dp)                :: Td(DTIMES(0:box%n_cell+1))
    real(dp)                :: dens((box%n_cell+2)**NDIM, n_species)
    real(dp)                :: rates((box%n_cell+2)**NDIM, n_reactions)
    real(dp)                :: derivs((box%n_cell+2)**NDIM, n_species)
    logical                 :: have_derivs

    if (.not. gas_constant_density) then
       N_inv = 1/box%cc(DTIMES(:), i_gas_dens)
    else
       N_inv = 1/gas_number_density
    end if

    Td = SI_to_Townsend * box%cc(DTIMES(:), i_electric_fld) * N_inv

    have_derivs = .false.
    nc = box%n_cell
    n_cells = (nc+2)**NDIM

    do n = 1, n_var
       select case (output_extra_vars(n))
       case ("eV")
          ! Add electron energy in eV
          new_vars(DTIMES(:), n) = LT_get_col(eV_vs_fld, 1, Td)
       case ("sigma")
          ! Add plasma conductivity (e mu n_e)
          new_vars(DTIMES(:), n) = LT_get_col(td_tbl, td_mobility, Td) * &
               N_inv * box%cc(DTIMES(:), i_electron) * UC_elem_charge
       case default
          ! Assume temporal production of some species is added, prefixed by "src_"

          if (.not. have_derivs) then
             call get_rates(pack(Td, .true.), rates, n_cells)

             dens(:, n_gas_species+1:n_species) = reshape(&
                  box%cc(DTIMES(:), species_itree(n_gas_species+1:n_species)), &
                  [n_cells, n_plasma_species])

             call get_derivatives(dens, rates, derivs, n_cells)
             have_derivs = .true.
          end if

          species_name = trim(output_extra_vars(n)(5:))
          i_species = species_index(species_name)

          if (i_species == -1) then
             print *, "No species corresponding to ", species_name
             error stop "Error adding time derivative to output"
          else
             new_vars(DTIMES(:), n) = reshape(derivs(:, i_species), &
                  [DTIMES(nc+2)])
          end if
       end select
    end do
  end subroutine add_variables

  subroutine output_log(tree, filename, out_cnt, wc_time)
    use m_field, only: current_voltage
    use m_user_methods
    use m_chemistry
    use m_analysis
    use m_dielectric
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt !< Output number
    real(dp), intent(in)         :: wc_time !< Wallclock time
    character(len=50), save      :: fmt
    integer                      :: my_unit, n
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(NDIM) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er, min_Er
    real(dp)                     :: sum_elem_charge, tmp, ne_zminmax(2)
    real(dp)                     :: elecdens_threshold, max_field_tip
    real(dp)                     :: r0(NDIM), r1(NDIM), r_tip(NDIM)
    type(af_loc_t)               :: loc_elec, loc_field, loc_Er, loc_tip
    integer                      :: i, n_reals, n_user_vars
    character(len=name_len)      :: var_names(user_max_log_vars)
    real(dp)                     :: var_values(user_max_log_vars)
    logical, save                :: first_time     = .true.

    if (associated(user_log_variables)) then
       ! Set default names for the user variables
       do i = 1, user_max_log_vars
          write(var_names, "(A,I0)") "uservar_", i
       end do
       call user_log_variables(tree, n_user_vars, var_names, var_values)
    else
       n_user_vars = 0
    end if

    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_1pos_ion, sum_pos_ion)
    call af_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    call af_tree_max_fc(tree, 1, electric_fld, max_Er, loc_Er)
    call af_tree_min_fc(tree, 1, electric_fld, min_Er)

    ! Scale threshold with gas number density
    elecdens_threshold = density_threshold * &
         (gas_number_density/2.414e25_dp)**2
    call analysis_zmin_zmax_threshold(tree, i_electron, elecdens_threshold, &
         [ST_domain_len(NDIM), 0.0_dp], ne_zminmax)

    ! Assume streamer tip is located at the farthest plasma density (above a
    ! threshold) away from the electrode boundaries
    r0 = ST_domain_origin
    r1 = ST_domain_origin + ST_domain_len

    if (ne_zminmax(1) - ST_domain_origin(NDIM) < &
         ST_domain_origin(NDIM) + ST_domain_len(NDIM) - ne_zminmax(2)) then
       r0(NDIM) = ne_zminmax(2) - 0.02_dp * ST_domain_len(NDIM)
       r1(NDIM) = ne_zminmax(2) + 0.02_dp * ST_domain_len(NDIM)
    else
       r0(NDIM) = ne_zminmax(1) - 0.02_dp * ST_domain_len(NDIM)
       r1(NDIM) = ne_zminmax(1) + 0.02_dp * ST_domain_len(NDIM)
    end if

    call analysis_max_var_region(tree, i_electric_fld, r0, r1, &
         max_field_tip, loc_tip)

    if (loc_tip%id > 0) then
       r_tip = af_r_loc(tree, loc_tip)
    else
       r_tip = 0.0_dp
    end if

    sum_elem_charge = 0
    do n = n_gas_species+1, n_species
       if (species_charge(n) /= 0) then
          call af_tree_sum_cc(tree, species_itree(n), tmp)
          sum_elem_charge = sum_elem_charge + tmp * species_charge(n)
       end if
    end do

    if (ST_use_dielectric) then
       call todo_diel_get_integral(diel, i_surf_dens, tmp)
       sum_elem_charge = sum_elem_charge + tmp
    end if

    dt = global_dt

    if (first_time) then
       first_time = .false.
       open(newunit=my_unit, file=trim(filename), action="write")
#if NDIM == 1
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &sum(charge) sum(J.E) max(E) x max(n_e) x voltage ne_zmin ne_zmax &
            &max(Etip) x wc_time n_cells min(dx) &
            &highest(lvl)"
#elif NDIM == 2
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &sum(charge) sum(J.E) max(E) x y max(n_e) x y max(E_r) x y min(E_r) voltage &
            &ne_zmin ne_zmax max(Etip) x y wc_time n_cells min(dx) highest(lvl)"
#elif NDIM == 3
       write(my_unit, "(A)", advance="no") "it time dt v sum(n_e) sum(n_i) &
            &sum(charge) sum(J.E) max(E) x y z max(n_e) x y z voltage &
            &ne_zmin ne_zmax max(Etip) x y z wc_time n_cells min(dx) highest(lvl)"
#endif
       if (associated(user_log_variables)) then
          do i = 1, n_user_vars
             write(my_unit, "(A)", advance="no") " "//trim(var_names(i))
          end do
       end if
       write(my_unit, *) ""
       close(my_unit)

       ! Start with velocity zero
       prev_pos = af_r_loc(tree, loc_field)
    end if

#if NDIM == 1
    n_reals = 17
#elif NDIM == 2
    n_reals = 24
#elif NDIM == 3
    n_reals = 23
#endif

    if (associated(user_log_variables)) then
       write(fmt, "(A,I0,A,I0,A)") "(I6,", n_reals, "E20.8,I12,1E20.8,I3,", &
            n_user_vars, "E20.8)"
    else
       write(fmt, "(A,I0,A)") "(I6,", n_reals, "E20.8,I12,1E20.8,I3)"
    end if

    velocity = norm2(af_r_loc(tree, loc_field) - prev_pos) / output_dt
    prev_pos = af_r_loc(tree, loc_field)

    open(newunit=my_unit, file=trim(filename), action="write", &
         position="append")
#if NDIM == 1
    write(my_unit, fmt) out_cnt, global_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_elem_charge, sum(ST_global_JdotE(1, :)), &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), current_voltage, ne_zminmax, &
         max_field_tip, r_tip, &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 2
    write(my_unit, fmt) out_cnt, global_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_elem_charge, sum(ST_global_JdotE(1, :)), &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), min_Er, &
         current_voltage, ne_zminmax, max_field_tip, r_tip, &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, global_time, dt, velocity, sum_elec, &
         sum_pos_ion, sum_elem_charge, sum(ST_global_JdotE(1, :)), &
         max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), current_voltage, ne_zminmax, &
         max_field_tip, r_tip, &
         wc_time, af_num_cells_used(tree), &
         af_min_dr(tree),tree%highest_lvl, &
         var_values(1:n_user_vars)
#endif
    close(my_unit)

  end subroutine output_log

  !> Write a net stoichiometric matrix to a file
  subroutine output_stoichiometric_matrix(filename)
    use m_chemistry
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, m, i, ix
    integer, allocatable         :: stoich(:,:)

    allocate(stoich(n_reactions, n_species))
    stoich(:, :) = 0

    do m = 1, n_reactions
       ! Subtract consumed species
       do i = 1, size(reactions(m)%ix_in)
          ix = reactions(m)%ix_in(i)
          stoich(m, ix) = stoich(m, ix) - 1
       end do

       ! Add produced species
       do i = 1, size(reactions(m)%ix_out)
          ix = reactions(m)%ix_out(i)
          stoich(m, ix) = stoich(m, ix) + &
               reactions(m)%multiplicity_out(i)
       end do
    end do

    open(newunit=my_unit, file=trim(filename), action="write")
    do i = 1, n_species
       write(my_unit, *) stoich(:, i)
    end do
    write(my_unit, *) ""
    close(my_unit)
  end subroutine output_stoichiometric_matrix

  !> Write a list of chemical species
  subroutine output_chemical_species(filename)
    use m_chemistry
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, i

    open(newunit=my_unit, file=trim(filename), action="write")
    do i = 1, n_species
       write(my_unit, "(A)") species_list(i)
    end do
    write(my_unit, "(A)") ""
    close(my_unit)
  end subroutine output_chemical_species

  !> Write a list of chemical reactions
  subroutine output_chemical_reactions(filename)
    use m_chemistry
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, i

    open(newunit=my_unit, file=trim(filename), action="write")
    do i = 1, n_reactions
       write(my_unit, "(A)") reactions(i)%description
    end do
    write(my_unit, "(A)") ""
    close(my_unit)
  end subroutine output_chemical_reactions

  !> Write space-integrated and time-integrated reaction rates
  subroutine output_chemical_rates(filename, first_time)
    use m_chemistry
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: first_time
    integer                      :: my_unit, iostate

    if (first_time) then
       ! Clear the file
       open(newunit=my_unit, file=trim(filename), iostat=iostate)
       if (iostate == 0) close(my_unit, status='delete')
    else
       open(newunit=my_unit, file=trim(filename), action="write", &
            position="append")
       write(my_unit, "(*(E20.8))") global_time, &
            sum(ST_global_rates(1:n_reactions, :), dim=2)
       close(my_unit)
    end if
  end subroutine output_chemical_rates

  !> Write space-integrated species densities
  subroutine output_chemical_amounts(tree, filename, first_time)
    use m_chemistry
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    logical, intent(in)          :: first_time
    integer                      :: my_unit, n, iostate
    real(dp)                     :: sum_dens(n_species)

    if (first_time) then
       ! Clear the file
       open(newunit=my_unit, file=trim(filename), iostat=iostate)
       if (iostate == 0) close(my_unit, status='delete')
    else
       do n = 1, n_species
          if (species_itree(n) > 0) then
             call af_tree_sum_cc(tree, species_itree(n), sum_dens(n))
          else
             sum_dens(n) = 0.0_dp ! A neutral gas specie
          end if
       end do

       open(newunit=my_unit, file=trim(filename), action="write", &
            position="append")
       write(my_unit, "(*(E20.8))") global_time, sum_dens
       close(my_unit)
    end if

  end subroutine output_chemical_amounts

  !> Write statistics to a file that can be used for regression testing
  subroutine output_regression_log(tree, filename, out_cnt, wc_time)
    use m_chemistry
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt !< Output number
    real(dp), intent(in)         :: wc_time !< Wallclock time
    character(len=30)            :: fmt
    integer                      :: my_unit, n
    real(dp)                     :: vol
    real(dp)                     :: sum_dens(n_species)
    real(dp)                     :: sum_dens_sq(n_species)
    real(dp)                     :: max_dens(n_species)

    vol = af_total_volume(tree)
    do n = 1, n_species
       if (species_itree(n) > 0) then
          call af_tree_sum_cc(tree, species_itree(n), sum_dens(n))
          call af_tree_sum_cc(tree, species_itree(n), sum_dens_sq(n), power=2)
          call af_tree_max_cc(tree, species_itree(n), max_dens(n))
       else
          ! A neutral gas specie
          sum_dens(n) = 0.0_dp
          sum_dens_sq(n) = 0.0_dp
          max_dens(n) = 0.0_dp
       end if
    end do

    if (out_cnt == 0) then
       open(newunit=my_unit, file=trim(filename), action="write")
       write(my_unit, "(A)", advance="no") "it time dt"
       do n = 1, n_species
          write(my_unit, "(A)", advance="no") &
               " sum(" // trim(species_list(n)) // ")"
       end do
       do n = 1, n_species
          write(my_unit, "(A)", advance="no") &
               " sum(" // trim(species_list(n)) // "^2)"
       end do
       do n = 1, n_species
          write(my_unit, "(A)", advance="no") &
               " max(" // trim(species_list(n)) // ")"
       end do
       write(my_unit, "(A)") ""
       close(my_unit)
    end if

    write(fmt, "(A,I0,A)") "(I0,", 3+3*n_species, "E20.8)"

    open(newunit=my_unit, file=trim(filename), action="write", &
         position="append")
    write(my_unit, fmt) out_cnt, global_time, global_dt, &
         sum_dens/vol, sum_dens_sq/vol, max_dens
    close(my_unit)

  end subroutine output_regression_log

  subroutine check_path_writable(pathname)
    character(len=*), intent(in) :: pathname
    integer                      :: my_unit, iostate
    open(newunit=my_unit, file=trim(pathname)//"_DUMMY", iostat=iostate)
    if (iostate /= 0) then
       print *, "Output name: " // trim(pathname) // '_...'
       error stop "Directory not writable (does it exist?)"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_path_writable

  !> Print short status message
  subroutine output_status(tree, time, wc_time, it, dt)
    use m_dt
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp), intent(in)   :: wc_time
    integer, intent(in)    :: it
    real(dp), intent(in)   :: dt
    write(*, "(F7.2,A,I0,A,E10.3,A,E10.3,A,E10.3,A,E10.3)") &
         100 * time / ST_end_time, "% it: ", it, &
         " t:", time, " dt:", dt, " wc:", wc_time, &
         " ncell:", real(af_num_cells_used(tree), dp)

    ! This line prints the different time step restrictions
    write(*, "(A,4E10.3,A)") "         dt: ", &
         minval(dt_matrix(1:dt_num_cond, :), dim=2), &
         " (cfl diff drt chem)"
  end subroutine output_status

  subroutine output_fld_maxima(tree, filename)
    use m_analysis
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer                      :: my_unit, i, n, n_found
    integer, parameter           :: n_max = 1000
    real(dp)                     :: coord_val(NDIM+1, n_max), d

    ! Get locations of the maxima
    call analysis_get_maxima(tree, i_electric_fld, field_maxima_threshold, &
         n_max, coord_val, n_found)

    if (n_found > n_max) then
       print *, "output_fld_maxima warning: n_found > n_max"
       print *, "n_found = ", n_found, "n_max = ", n_max
       n_found = n_max
    end if

    ! Merge maxima that are too close
    do n = n_found, 1, -1
       do i = 1, n-1
          d = norm2(coord_val(1:NDIM, i) - coord_val(1:NDIM, n))
          if (d < field_maxima_distance) then
             ! Replace the smaller value by the one at the end of the list, and
             ! shorten the list by one item
             if (coord_val(NDIM+1, i) < coord_val(NDIM+1, n)) then
                coord_val(:, i) = coord_val(:, n)
             end if
             coord_val(:, n) = coord_val(:, n_found)
             n_found = n_found - 1
             exit
          end if
       end do
    end do

    open(newunit=my_unit, file=trim(filename), action="write")
    do n = 1, n_found
       if (coord_val(NDIM+1, n) > field_maxima_threshold) then
          write(my_unit, *) coord_val(:, n)
       end if
    end do
    close(my_unit)

  end subroutine output_fld_maxima

#if NDIM == 2
  subroutine output_cross(tree, filename)
    use m_af_all
    use m_analysis
    use m_streamer
    type(af_t), intent(in) :: tree
    character(len=*), intent(inout) :: filename

    integer  :: my_unit
    integer  :: i
    real(dp) :: z, elec_dens, charge_dens, current_dens

    open(newunit=my_unit, file=trim(filename)//".txt", action="write")
    write(my_unit, '(A)') "z elec_dens charge_dens current_dens"

    do i = 1, cross_npoints
       z = i * ST_domain_len(2) / (cross_npoints + 1)
       call analysis_get_cross(tree, cross_rmax, z, elec_dens, &
            charge_dens, current_dens)
       write(my_unit, *) z, elec_dens, charge_dens, current_dens
    end do

    close(my_unit)
  end subroutine output_cross
#endif

  subroutine set_power_density_box(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: J_dot_E

    nc = box%n_cell
    do KJI_DO(1, nc)
       ! Compute inner product flux * field over the cell faces
       J_dot_E = 0.5_dp * sum(box%fc(IJK, :, flux_elec) * box%fc(IJK, :, electric_fld))
#if NDIM == 1
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, 1, flux_elec) * box%fc(i+1, 1, electric_fld))
#elif NDIM == 2
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, j, 1, flux_elec) * box%fc(i+1, j, 1, electric_fld) + &
            box%fc(i, j+1, 2, flux_elec) * box%fc(i, j+1, 2, electric_fld))
#elif NDIM == 3
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, j, k, 1, flux_elec) * box%fc(i+1, j, k, 1, electric_fld) + &
            box%fc(i, j+1, k, 2, flux_elec) * box%fc(i, j+1, k, 2, electric_fld) + &
            box%fc(i, j, k+1, 3, flux_elec) * box%fc(i, j, k+1, 3, electric_fld))
#endif

       box%cc(IJK, i_power_density) = J_dot_E * UC_elec_charge
    end do; CLOSE_DO
  end subroutine set_power_density_box

  subroutine set_gas_primitives_box(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc, idim

    nc = box%n_cell
    do KJI_DO(1, nc)
       do idim = 1, NDIM
          ! Compute velocity components
          box%cc(IJK, gas_prim_vars(i_mom(idim))) = &
               box%cc(IJK, gas_vars(i_mom(idim))) / box%cc(IJK, gas_vars(i_rho))
       end do

       ! Compute the pressure
       box%cc(IJK, gas_prim_vars(i_e)) = (gas_euler_gamma-1.0_dp) * (box%cc(IJK,gas_vars( i_e)) - &
            0.5_dp*box%cc(IJK,gas_vars( i_rho))* sum(box%cc(IJK, gas_vars(i_mom(:)))**2))
       ! Compute the temperature (T = P/(n*kB), n=gas number density)
       box%cc(IJK, gas_prim_vars(i_e+1)) = box%cc(IJK, gas_prim_vars(i_e))/ &
            (box%cc(IJK, i_gas_dens)* UC_boltzmann_const)

    end do; CLOSE_DO
  end subroutine set_gas_primitives_box

  !> @todo replace this by routine in afivo/src/m_dielectric
  subroutine todo_diel_get_integral(diel, i_surf, surf_int)
    type(surfaces_t), intent(inout) :: diel
    integer, intent(in)               :: i_surf   !< Surface variables
    real(dp), intent(out)             :: surf_int !< Surface integral
    integer                           :: ix

    surf_int = 0
    do ix = 1, diel%max_ix
       if (diel%surfaces(ix)%in_use) then
#if NDIM == 2
          surf_int = surf_int + product(diel%surfaces(ix)%dr) * &
               sum(diel%surfaces(ix)%sd(:, i_surf))
#elif NDIM == 3
          surf_int = surf_int + product(diel%surfaces(ix)%dr) * &
               sum(diel%surfaces(ix)%sd(:, :, i_surf))
#endif
       end if
    end do
  end subroutine todo_diel_get_integral

end module m_output
