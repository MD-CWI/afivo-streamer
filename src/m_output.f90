!> Module with methods and parameters related to writing output
module m_output
  use m_af_all
  use m_types
  use m_advance, only: advance_max_dt
  use m_streamer

  implicit none
  private

  ! Name for the output files
  character(len=string_len), public, protected :: output_name = "output/my_sim"

  ! Optional variable (to show ionization source term)
  integer, public, protected :: i_src = -1 ! Source term

  ! If true, only include n_e, n_i and |E| in output files
  logical, public, protected :: output_compact = .false.

  ! Time between writing output
  real(dp), public, protected :: output_dt = 1.0e-10_dp

  ! Include ionization source term in output
  logical, public, protected :: output_src_term = .false.

  ! If positive: decay rate for source term (1/s) for time-averaged values
  real(dp), public, protected :: output_src_decay_rate = -1.0_dp

  ! Write output along a line
  logical, public, protected :: lineout_write = .false.

  ! Write Silo output
  logical, public, protected :: silo_write = .true.

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

  ! Write uniform output in a plane
  logical, public, protected :: plane_write = .false.

  ! Which variable to include in plane
  integer, public, protected :: plane_ivar

  ! Use this many points for plane data
  integer, public, protected :: plane_npixels(2) = [64, 64]

  ! Relative position of plane minimum coordinate
  real(dp), public, protected :: plane_rmin(NDIM) = 0.0_dp

  ! Relative position of plane maximum coordinate
  real(dp), public, protected :: plane_rmax(NDIM) = 1.0_dp

  ! Print status every this many seconds
  real(dp), public, protected :: output_status_delay = 60.0_dp

  ! Public methods
  public :: output_initialize
  public :: output_write
  public :: output_log
  public :: output_status

contains

  subroutine output_initialize(tree, cfg)
    use m_config
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    character(len=name_len)    :: varname
    real(dp)                   :: tmp

    call CFG_add_get(cfg, "output%name", output_name, &
         "Name for the output files (e.g. output/my_sim)")

    call check_path_writable(output_name)

    call CFG_add_get(cfg, "output%status_delay", output_status_delay, &
         "Print status every this many seconds")

    call CFG_add_get(cfg, "output%compact", output_compact, &
         "If true, only include n_e, n_i and |E| in output files")

    if (output_compact) then
       tree%cc_write_output(:) = .false.
       tree%cc_write_output([i_electron, i_1pos_ion, i_electric_fld]) = .true.
    end if

    call CFG_add_get(cfg, "output%dt", output_dt, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "output%src_term", output_src_term, &
         "Include ionization source term in output")

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

    call CFG_add_get(cfg, "datfile%write", datfile_write, &
         "Write binary output files (to resume later)")
    call CFG_add_get(cfg, "datfile%per_outputs", datfile_per_outputs, &
         "Write binary output files every N outputs")

    call CFG_add_get(cfg, "lineout%write", lineout_write, &
         "Write output along a line")
    call CFG_add_get(cfg, "lineout%npoints", lineout_npoints, &
         "Use this many points for lineout data")
    call CFG_add_get(cfg, "lineout%rmin", lineout_rmin(1:NDIM), &
         "Relative position of line minimum coordinate")
    call CFG_add_get(cfg, "lineout%rmax", lineout_rmax(1:NDIM), &
         "Relative position of line maximum coordinate")

    call CFG_add_get(cfg, "plane%write", plane_write, &
         "Write uniform output in a plane")
    varname = "e"
    call CFG_add_get(cfg, "plane%varname", varname, &
         "Names of variable to write in a plane")
    plane_ivar = af_find_cc_variable(tree, trim(varname))
    call CFG_add_get(cfg, "plane%npixels", plane_npixels, &
         "Use this many pixels for plane data")
    call CFG_add_get(cfg, "plane%rmin", plane_rmin(1:NDIM), &
         "Relative position of plane minimum coordinate")
    call CFG_add_get(cfg, "plane%rmax", plane_rmax(1:NDIM), &
         "Relative position of plane maximum coordinate")

  end subroutine output_initialize

  subroutine output_write(tree, output_cnt, wc_time, write_sim_data)
    use m_photoi
    use m_field
    use m_user_methods
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

    if (silo_write) then
       ! Because the mesh could have changed
       if (photoi_enabled) call photoi_set_src(tree, advance_max_dt)
       call field_set_rhs(tree, 0)

       do i = 1, tree%n_var_cell
          if (tree%cc_write_output(i)) then
             call af_restrict_tree(tree, i)
             call af_gc_tree(tree, i)
          end if
       end do

       write(fname, "(A,I6.6)") trim(output_name) // "_", output_cnt
       call af_write_silo(tree, fname, output_cnt, global_time)
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

    if (plane_write) then
       write(fname, "(A,I6.6)") trim(output_name) // &
            "_plane_", output_cnt
       call af_write_plane(tree, fname, [plane_ivar], &
            plane_rmin * ST_domain_len + ST_domain_origin, &
            plane_rmax * ST_domain_len + ST_domain_origin, &
            plane_npixels)
    end if

    if (lineout_write) then
       write(fname, "(A,I6.6)") trim(output_name) // &
            "_line_", output_cnt
       call af_write_line(tree, trim(fname), &
            [i_electron, i_1pos_ion, i_phi, i_electric_fld], &
            r_min = lineout_rmin(1:NDIM) * ST_domain_len + ST_domain_origin, &
            r_max = lineout_rmax(1:NDIM) * ST_domain_len + ST_domain_origin, &
            n_points=lineout_npoints)
    end if

  end subroutine output_write

  subroutine output_log(tree, filename, out_cnt, wc_time)
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt
    real(dp), intent(in)         :: wc_time
    character(len=30), save      :: fmt
    integer, parameter           :: my_unit        = 123
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(NDIM) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er, min_Er
    type(af_loc_t)               :: loc_elec, loc_field, loc_Er

    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_1pos_ion, sum_pos_ion)
    call af_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    call af_tree_max_fc(tree, 1, electric_fld, max_Er, loc_Er)
    call af_tree_min_fc(tree, 1, electric_fld, min_Er)

    dt = advance_max_dt

    if (out_cnt == 1) then
       open(my_unit, file=trim(filename), action="write")
#if NDIM == 2
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y max(n_e) x y max(E_r) x y min(E_r) wc_time n_cells min(dx) highest(lvl)"
#elif NDIM == 3
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y z max(n_e) x y z wc_time n_cells min(dx) highest(lvl)"
#endif
       close(my_unit)

       ! Start with velocity zero
       prev_pos = af_r_loc(tree, loc_field)
    end if

#if NDIM == 2
    fmt = "(I6,16E16.8,I12,1E16.8,I3)"
#elif NDIM == 3
    fmt = "(I6,14E16.8,I12,1E16.8,I3)"
#endif

    velocity = norm2(af_r_loc(tree, loc_field) - prev_pos) / output_dt
    prev_pos = af_r_loc(tree, loc_field)

    open(my_unit, file=trim(filename), action="write", &
         position="append")
#if NDIM == 2
    write(my_unit, fmt) out_cnt, global_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), min_Er, &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, global_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl
#endif
    close(my_unit)

  end subroutine output_log

  subroutine check_path_writable(pathname)
    character(len=*), intent(in) :: pathname
    integer                      :: my_unit, iostate
    open(newunit=my_unit, file=trim(pathname)//"_DUMMY", iostat=iostate)
    if (iostate /= 0) then
       print *, "Output directory: " // trim(pathname)
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

end module m_output
