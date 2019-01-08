#include "../afivo/src/cpp_macros.h"
!> Module to compute electric fields
module m_field
  use m_af_all
  use m_types
  use m_streamer

  implicit none
  private

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
  real(dp) :: field_voltage

  !> Drop-off radius
  real(dp) :: field_dropoff_radius = 1e-3_dp

  !> Relative width over which the potential drops
  real(dp) :: field_dropoff_relwidth = 0.5_dp

  !> Location from which the field drops off (set below)
  real(dp) :: field_dropoff_pos(2) = 0.0_dp

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
  public :: field_from_potential
  public :: field_get_amplitude
  public :: field_set_voltage

  public :: field_bc_homogeneous
  public :: field_bc_dropoff_lin
  public :: field_bc_dropoff_log

contains

  !> Initialize this module
  subroutine field_initialize(tree, cfg, mg)
    use m_config
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg !< Settings
    type(mg_t), intent(inout)  :: mg  !< Multigrid option struct

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
    case ("point_charge")
       mg%sides_bc => field_bc_point_charge
    case default
       error stop "field_bc_select error: invalid condition"
    end select

    ! Set the multigrid options. First define the variables to use
    mg%i_phi = i_phi
    mg%i_tmp = i_electric_fld
    mg%i_rhs = i_rhs
    if (ST_use_dielectric) mg%i_eps = i_eps

    ! This automatically handles cylindrical symmetry
    mg%box_op => mg_auto_op
    mg%box_gsrb => mg_auto_gsrb
    mg%box_corr => mg_auto_corr

    ! This routine always needs to be called when using multigrid
    call mg_init_mg(mg)

    call af_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)
    call af_set_cc_methods(tree, i_electric_fld, &
         af_bc_neumann_zero, af_gc_interp)

  end subroutine field_initialize

  !> Compute electric field on the tree. First perform multigrid to get electric
  !> potential, then take numerical gradient to geld field.
  subroutine field_compute(tree, mg, s_in, time, have_guess)
    use m_units_constants
    use m_chemistry
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg ! Multigrid option struct
    integer, intent(in)       :: s_in
    real(dp), intent(in)      :: time
    logical, intent(in)       :: have_guess
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

    call field_set_voltage(tree, time)

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
  function field_get_amplitude(tree, time) result(electric_fld)
    use m_units_constants
    type(af_t), intent(in) :: tree
    real(dp), intent(in)    :: time
    real(dp)                :: electric_fld, t_rel
    type(af_loc_t)         :: loc_field
    real(dp)                :: r(NDIM), zrel, max_fld

    t_rel = time - field_mod_t0
    t_rel = min(t_rel, field_mod_t1-field_mod_t0)

    if (t_rel > 0) then
       electric_fld = field_amplitude * exp(-t_rel/field_decay_time) + &
            t_rel * field_lin_deriv + &
            field_sin_amplitude * &
            sin(t_rel * field_sin_freq * 2 * UC_pi)
    else if (field_stability_search) then
       call af_tree_max_cc(tree, i_electric_fld, max_fld, loc_field)
       r = af_r_loc(tree, loc_field)
       zrel = r(NDIM) / ST_domain_len
       zrel = (zrel - field_stability_zmin) / &
            (field_stability_zmax - field_stability_zmin)

       if (zrel > 0.0_dp .and. max_fld > field_stability_threshold) then
          electric_fld = (1 - zrel) * field_amplitude
       else
          electric_fld = field_amplitude
       end if
    else
       electric_fld = field_amplitude
    end if

  end function field_get_amplitude

  !> Compute the voltage at a given time
  subroutine field_set_voltage(tree, time)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)    :: time

    field_voltage = -ST_domain_len * field_get_amplitude(tree, time)
  end subroutine field_set_voltage

  !> This fills ghost cells near physical boundaries for the potential
  subroutine field_bc_homogeneous(box, nb, iv, bc_type)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, iv) = 0
    case (af_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0
    case (af_neighb_lowy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc,    0, iv) = 0
    case (af_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = field_voltage
#elif NDIM == 3
    case (af_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, 1:nc, iv) = 0
    case (af_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, 1:nc, iv) = 0
    case (af_neighb_lowy)
       bc_type = af_bc_neumann
       box%cc(1:nc,    0, 1:nc, iv) = 0
    case (af_neighb_highy)
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, 1:nc, iv) = 0
    case (af_neighb_lowz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc,    0, iv) = 0
    case (af_neighb_highz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc, nc+1, iv) = field_voltage
#endif
    end select

  end subroutine field_bc_homogeneous

  subroutine field_bc_dropoff_lin(box, nb, iv, bc_type)
    type(box_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc, i
#if NDIM == 3
    integer                      :: j
#endif
    real(dp)                     :: rr(NDIM), rdist

    nc = box%n_cell

    select case (nb)
#if NDIM == 2
    case (af_neighb_highy)
       bc_type = af_bc_dirichlet

       do i = 1, nc
          rr = af_r_cc(box, [i, 0])
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
#elif NDIM == 3
    case (af_neighb_highz)
       bc_type = af_bc_dirichlet

       do j = 1, nc
          do i = 1, nc
             rr = af_r_cc(box, [i, j, 0])
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
    type(box_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc, i
#if NDIM == 3
    integer                      :: j
#endif
    real(dp)                     :: rr(NDIM), rdist, tmp

    nc = box%n_cell
    tmp = field_dropoff_relwidth * ST_domain_len

    select case (nb)
#if NDIM == 2
    case (af_neighb_highy)
       bc_type = af_bc_dirichlet

       do i = 1, nc
          rr = af_r_cc(box, [i, 0])
          rdist = abs(rr(1) - field_dropoff_pos(1))

          if (rdist < field_dropoff_radius) then
             box%cc(i, nc+1, iv) = field_voltage
          else
             box%cc(i, nc+1, iv) = field_voltage * &
                  log(1 + tmp/rdist) / log(1 + tmp/field_dropoff_radius)
          end if
       end do
#elif NDIM == 3
    case (af_neighb_highz)
       bc_type = af_bc_dirichlet

       do j = 1, nc
          do i = 1, nc
             rr = af_r_cc(box, [i, j, 0])
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

  !> Create a field of the form E = E_0 - c / r^2
  subroutine field_bc_point_charge(box, nb, iv, bc_type)
    use m_units_constants
    type(box_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc, i, j
    real(dp)                     :: rr(NDIM), r0(NDIM), r1(NDIM), q

    nc = box%n_cell
    bc_type = af_bc_dirichlet
    q = field_point_charge / (4 * UC_pi * UC_eps0)
    r0 = field_point_r0 * ST_domain_len
    r1 = r0
    r1(NDIM) = -r0(NDIM)

    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
       if (ST_cylindrical) then
          bc_type = af_bc_neumann
          box%cc(0, 1:nc, iv) = 0
       else
          do j = 1, nc
             rr = af_rr_cc(box, [0.5_dp, j+0.0_dp])
             box%cc(0, j, iv) = q/norm2(rr-r0) - q/norm2(rr-r1)
          end do
       end if
    case (af_neighb_highx)
       do j = 1, nc
          rr = af_rr_cc(box, [nc+0.5_dp, j+0.0_dp])
          box%cc(nc+1, j, iv) = q/norm2(rr-r0) - q/norm2(rr-r1)
       end do
    case (af_neighb_lowy)
       do i = 1, nc
          rr = af_rr_cc(box, [i+0.0_dp, 0.5_dp])
          box%cc(i, 0, iv) = q/norm2(rr-r0) - q/norm2(rr-r1)
       end do
    case (af_neighb_highy)
       do i = 1, nc
          rr = af_rr_cc(box, [i+0.0_dp, nc+0.5_dp])
          box%cc(i, nc+1, iv) = q/norm2(rr-r0) - q/norm2(rr-r1)
       end do
#elif NDIM == 3
#endif
    case default
       error stop "Not implemented"
       i = 0; j = 0; rr = 0;
       call field_bc_homogeneous(box, nb, iv, bc_type)
    end select
  end subroutine field_bc_point_charge

  !> Compute electric field from electrical potential
  subroutine field_from_potential(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    real(dp)                   :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if NDIM == 2
    box%fc(1:nc+1, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi)) + &
         field_background(1)
    box%fc(1:nc, 1:nc+1, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi)) + &
         field_background(2)

    if (ST_use_dielectric) then
       ! Compute fields at the boundaries of the box, where eps can change
       box%fc(1, 1:nc, 1, electric_fld) = 2 * inv_dr * &
            (box%cc(0, 1:nc, i_phi) - box%cc(1, 1:nc, i_phi)) * &
            box%cc(0, 1:nc, i_eps) / &
            (box%cc(1, 1:nc, i_eps) + box%cc(0, 1:nc, i_eps))
       box%fc(nc+1, 1:nc, 1, electric_fld) = 2 * inv_dr * &
            (box%cc(nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, i_phi)) * &
            box%cc(nc+1, 1:nc, i_eps) / &
            (box%cc(nc+1, 1:nc, i_eps) + box%cc(nc, 1:nc, i_eps))
       box%fc(1:nc, 1, 2, electric_fld) = 2 * inv_dr * &
            (box%cc(1:nc, 0, i_phi) - box%cc(1:nc, 1, i_phi)) * &
            box%cc(1:nc, 0, i_eps) / &
            (box%cc(1:nc, 1, i_eps) + box%cc(1:nc, 0, i_eps))
       box%fc(1:nc, nc+1, 2, electric_fld) = 2 * inv_dr * &
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
    box%fc(1:nc+1, 1:nc, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, 1:nc, i_phi) - &
         box%cc(1:nc+1, 1:nc, 1:nc, i_phi)) + &
         field_background(1)
    box%fc(1:nc, 1:nc+1, 1:nc, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, 1:nc, i_phi) - &
         box%cc(1:nc, 1:nc+1, 1:nc, i_phi)) + &
         field_background(2)
    box%fc(1:nc, 1:nc, 1:nc+1, 3, electric_fld) = inv_dr * &
         (box%cc(1:nc, 1:nc, 0:nc, i_phi) - &
         box%cc(1:nc, 1:nc, 1:nc+1, i_phi)) + &
         field_background(3)

    if (ST_use_dielectric) then
       ! Compute fields at the boundaries of the box, where eps can change
       box%fc(1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr * &
            (box%cc(0, 1:nc, 1:nc, i_phi) - box%cc(1, 1:nc, 1:nc, i_phi)) * &
            box%cc(0, 1:nc, 1:nc, i_eps) / &
            (box%cc(1, 1:nc, 1:nc, i_eps) + box%cc(0, 1:nc, 1:nc, i_eps))
       box%fc(nc+1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr * &
            (box%cc(nc, 1:nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, 1:nc, i_phi)) * &
            box%cc(nc+1, 1:nc, 1:nc, i_eps) / &
            (box%cc(nc+1, 1:nc, 1:nc, i_eps) + box%cc(nc, 1:nc, 1:nc, i_eps))
       box%fc(1:nc, 1, 1:nc, 2, electric_fld) = 2 * inv_dr * &
            (box%cc(1:nc, 0, 1:nc, i_phi) - box%cc(1:nc, 1, 1:nc, i_phi)) * &
            box%cc(1:nc, 0, 1:nc, i_eps) / &
            (box%cc(1:nc, 1, 1:nc, i_eps) + box%cc(1:nc, 0, 1:nc, i_eps))
       box%fc(1:nc, nc+1, 1:nc, 2, electric_fld) = 2 * inv_dr * &
            (box%cc(1:nc, nc, 1:nc, i_phi) - box%cc(1:nc, nc+1, 1:nc, i_phi)) * &
            box%cc(1:nc, nc+1, 1:nc, i_eps) / &
            (box%cc(1:nc, nc+1, 1:nc, i_eps) + box%cc(1:nc, nc, 1:nc, i_eps))
       box%fc(1:nc, 1:nc, 1, 3, electric_fld) = 2 * inv_dr * &
            (box%cc(1:nc, 1:nc, 0, i_phi) - box%cc(1:nc, 1:nc, 1, i_phi)) * &
            box%cc(1:nc, 1:nc, 0, i_eps) / &
            (box%cc(1:nc, 1:nc, 1, i_eps) + box%cc(1:nc, 1:nc, 0, i_eps))
       box%fc(1:nc, 1:nc, nc+1, 3, electric_fld) = 2 * inv_dr * &
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
