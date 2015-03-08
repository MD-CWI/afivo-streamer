!> \example This example shows how one can do a more complicated simulation of
!> streamer discharges. This involves advection, diffusion and reaction
!> coupled to a Poisson equation. Furthermore, a dielectric material can be
!> placed in the domain.

program test_streamer_2d

  ! We need to import afivo, multigrid, and extra routines for dielectrics
  use m_afivo_2d
  use m_mg_2d
  use m_write_silo

  implicit none

  integer, parameter :: dp = kind(0.0d0)

  ! The size of the boxes that we use to construct our mesh
  integer, parameter :: box_size   = 8

  ! Box types, these are used to distinguish between different regions
  integer, parameter :: full_diel_box = 1
  integer, parameter :: part_diel_box = 2
  integer, parameter :: full_gas_box = 3

  ! Indices of cell-centered variables
  integer, parameter :: n_var_cell = 9
  integer, parameter :: i_elec     = 1 ! Electron density
  integer, parameter :: i_pion     = 2 ! Positive ion density
  integer, parameter :: i_elec_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pion_old = 4 ! For time-stepping scheme
  integer, parameter :: i_phi      = 5 ! Electrical potential
  integer, parameter :: i_fld      = 6 ! Electric field norm
  integer, parameter :: i_rhs      = 7 ! Source term Poisson
  integer, parameter :: i_res      = 8 ! Residual (multigrid)
  integer, parameter :: i_eps      = 9 ! Dielectric permittivity
  character(len=10)  :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs", "res", "eps"]

  ! Indices of face-centered variables
  integer, parameter :: n_var_face = 2
  integer, parameter :: f_elec     = 1 ! Electron flux
  integer, parameter :: f_fld      = 2 ! Electric field vector

  type(a2_t)         :: tree    ! This contains the full grid information
  type(mg2_t)        :: mg      ! Multigrid option struct
  type(ref_info_t)   :: ref_info

  integer            :: i, n, n_steps
  integer            :: output_cnt
  real(dp)           :: dt, time, end_time
  real(dp)           :: dt_adapt, dt_output
  character(len=40)  :: fname
  logical            :: write_out

  ! How many multigrid FMG cycles we perform per time step
  integer, parameter :: n_fmg_cycles = 1

  ! The (constant) diffusion coefficient for the electron density
  real(dp), parameter :: diff_coeff = 0.1_dp ! m^2/s

  ! The (constant) mobility of the electron density
  real(dp), parameter :: mobility = -0.03_dp ! m^2/(Vs)

  ! The length of the domain in each direction
  real(dp), parameter :: domain_len = 32.0e-3_dp ! m

  ! The location of the initial seed
  real(dp), parameter :: seed_r0(2) = &
       [0.35_dp, 0.5_dp] * domain_len - [0, 1] * 1e-3_dp
  real(dp), parameter :: seed_r1(2) = &
       [0.35_dp, 0.5_dp] * domain_len + [0, 1] * 1e-3_dp

  ! We store some extra data per box (its type)
  type :: my_data
     integer :: box_type
  end type my_data

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs
  mg%i_res        = i_res
  mg%i_eps        = i_eps

  ! The number of cycles at the lowest level
  mg%n_cycle_base = 8

  ! Routines to use for ...
  mg%sides_bc     => sides_bc_pot ! Filling ghost cell on physical boundaries
  mg%box_op       => auto_box_op  ! Performing the Laplacian
  mg%box_gsrb     => auto_box_gsrb ! Performing Gauss-Seidel relaxation
  mg%box_corr     => auto_box_corr ! Correcting the children of a grid

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  ! Set up the initial conditions
  do i = 1, 10
     call a2_loop_box(tree, set_init_cond)
     call a2_loop_box(tree, set_box_type)
     call a2_restrict_tree(tree, i_rhs)
     call compute_fld(tree, n_fmg_cycles)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  call a2_loop_box(tree, set_init_cond)

  output_cnt       = 0          ! Number of output files written
  time             = 0          ! Simulation time (all times are in s)
  dt_adapt         = 1.0e-11_dp ! Time per adaptation of the grid
  dt_output        = 2.5e-10_dp ! Time between writing output
  end_time         = 20.0e-9_dp ! Time to stop the simulation

  do
     ! Get a new time step, which is at most dt_adapt
     dt      = get_max_dt(tree)
     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?"
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_cnt * dt_output < time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I0,A)") "test_str2d_", output_cnt, ".silo"
     else
        write_out = .false.
     end if

     ! if (write_out) call a2_write_vtk(tree, trim(fname), &
     !      cc_names, output_cnt, time)
     if (write_out) call a2_write_silo(tree, trim(fname), &
          cc_names, output_cnt, time)

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, n_steps
        time = time + dt

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! First calculate fluxes
           call a2_loop_boxes(tree, fluxes_koren, .true.)

           call compute_fld(tree, n_fmg_cycles)

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_elec)
           call a2_restrict_tree(tree, i_pion)

           ! Fill ghost cells
           call a2_gc_sides(tree, i_elec, a2_sides_interp, sides_bc_dens)
           call a2_gc_sides(tree, i_pion, a2_sides_interp, sides_bc_dens)
        end do

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_dens, .true.)
     end do

     call a2_adjust_refinement(tree, set_ref_flags, ref_info)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_children(tree, ref_info)

        ! Compute the field on the
        call compute_fld(tree, n_fmg_cycles)

        ! This will every now-and-then clean up the data in the tree
        call a2_tidy_up(tree, 0.9_dp, 0.5_dp, 5000, .false.)
     end if
  end do

  call a2_destroy(tree)

contains

  ! Initialize the AMR tree
  subroutine init_tree(tree)
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 10*1000

    dr = domain_len / box_size

    ! Initialize tree
    call a2_init(tree, box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=4, n_boxes = n_boxes_init)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...
    nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

    ! Create the base mesh
    call a2_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! Refinement function
  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: nc
    real(dp)                 :: crv_phi, dr2, max_edens, max_fld

    nc        = boxes(id)%n_cell
    dr2       = boxes(id)%dr**2
    crv_phi   = dr2 * maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs)))
    max_edens = maxval(boxes(id)%cc(1:nc, 1:nc, i_elec))
    max_fld   = maxval(boxes(id)%cc(1:nc, 1:nc, i_fld))

    if (boxes(id)%dr > 2e-3_dp .or. (boxes(id)%dr > 1e-4_dp .and. &
         (a2_r_inside(boxes(id), seed_r0, 1.0e-3_dp) .or. &
         a2_r_inside(boxes(id), seed_r1, 1.0e-3_dp)))) then
       ref_flags(id) = a5_do_ref
    else if (crv_phi > 2.0e1_dp .and. max_fld > 3e6_dp &
         .and. boxes(id)%dr > 5e-6_dp) then
       ref_flags(id) = a5_do_ref
    else if (crv_phi < 2.0_dp) then
       ref_flags(id) = a5_rm_ref
    else
       ref_flags(id) = a5_kp_ref
    end if
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), sigma
    real(dp)                    :: seed_dens, bg_dens, dens

    nc        = box%n_cell
    ! The initial electron/ion seed lies between xy0 and xy1
    sigma     = 2e-4_dp          ! Width of the inital seed
    bg_dens   = 5e13             ! Background density
    seed_dens = 5e19

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])
          dens = seed_dens * rod_dens(xy, seed_r0, seed_r1, sigma, 3)

          if (xy(2) < 0.25_dp * domain_len) then
             box%cc(i, j, i_eps) = 5.0_dp
             box%cc(i, j, i_elec) = 0
             box%cc(i, j, i_pion) = 0
          else
             box%cc(i, j, i_eps) = 1.0_dp
             box%cc(i, j, i_elec) = bg_dens + dens
             box%cc(i, j, i_pion) = bg_dens + dens
          end if

       end do
    end do

    box%cc(:, :, i_phi) = 0     ! Inital potential set to zero
  end subroutine set_init_cond

  ! Store my_data (which contains the box type) for a box
  subroutine set_box_type(box)
    type(box2_t), intent(inout) :: box
    real(dp)                    :: max_eps, min_eps
    type(my_data)               :: md

    if (allocated(box%ud)) return ! Already set

    max_eps = maxval(box%cc(:,:, i_eps))
    min_eps = minval(box%cc(:,:, i_eps))

    if (max_eps > 1.0_dp) then
       if (min_eps > 1.0_dp) then
          md%box_type = full_diel_box
       else
          md%box_type = part_diel_box
       end if
    else
       md%box_type = full_gas_box
    end if

    box%ud = transfer(md, box%ud)
  end subroutine set_box_type

  ! Get maximum time step based on e.g. CFL criteria
  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    real(dp)               :: max_fld, max_dns, dr_min
    real(dp)               :: dt_cfl, dt_dif, dt_drt, dt_alpha

    call a2_tree_max_cc(tree, i_fld, max_fld)
    call a2_tree_max_cc(tree, i_elec, max_dns)

    dr_min = a2_min_dr(tree)

    ! CFL condition
    dt_cfl = dr_min / (abs(mobility) * max_fld)

    ! Diffusion condition
    dt_dif = dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * abs(mobility) * max_dns)

    ! Ionization limit
    dt_alpha =  1 / (abs(mobility) * max_fld * &
         max(epsilon(1.0_dp), get_alpha(max_fld)))

    get_max_dt = 0.75_dp * min(1/(1/dt_cfl + 1/dt_dif), dt_drt, dt_alpha)
  end function get_max_dt

  ! Used to create initial electron/ion seed
  real(dp) function rod_dens(xy, xy0, xy1, sigma, falloff_type)
    real(dp), intent(in) :: xy(2), xy0(2), xy1(2), sigma
    integer, intent(in)  :: falloff_type
    real(dp)             :: distance, temp

    distance = dist_line(xy, xy0, xy1)

    select case (falloff_type)
    case (1)                    ! Sigmoid
       rod_dens    = 2 / (1 + exp(distance / sigma))
    case (2)                    ! Gaussian
       rod_dens    = exp(-(distance/sigma)**2)
    case (3)                    ! Smooth-step
       if (distance < sigma) then
          rod_dens = 1
       else if (distance < 2 * sigma) then
          temp = distance/sigma - 1
          rod_dens = (1- (3 * temp**2 - 2 * temp**3))
       else
          rod_dens = 0.0_dp
       end if
    case default
       rod_dens = 0.0_dp
    end select
  end function rod_dens

  real(dp) function dist_line(xy, xy0, xy1)
    real(dp), intent(in) :: xy(2), xy0(2), xy1(2)
    real(dp) :: line_len2, temp
    real(dp) :: projection(2)

    line_len2 = sum((xy1 - xy0)**2)
    temp = sum((xy - xy0) * (xy1 - xy0)) / line_len2

    if (temp < 0.0_dp) then
       dist_line = sqrt(sum((xy-xy0)**2))
    else if (temp > 1.0_dp) then
       dist_line = sqrt(sum((xy-xy1)**2))
    else
       projection = xy0 + temp * (xy1 - xy0)
       dist_line = sqrt(sum((xy-projection)**2))
    end if
  end function dist_line

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, n_fmg)
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: n_fmg

    real(dp), parameter :: UC_eps0 = 8.8541878176d-12
    real(dp), parameter :: UC_elem_charge = 1.6022d-19
    real(dp), parameter :: fac = UC_elem_charge / UC_eps0
    integer :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(:, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, i_elec) - &
               tree%boxes(id)%cc(:, :, i_pion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Restrict the rhs
    call a2_restrict_tree(tree, i_rhs)

    ! Perform n_fmg full-multigrid cycles
    do i = 1, n_fmg
       call mg2_fas_fmg(tree, mg)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)
  end subroutine compute_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, f_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, f_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    ! Compute fields at the boundaries of the box, where eps can change
    box%fx(1, :, f_fld) = 2 * inv_dr * &
         (box%cc(0, 1:nc, i_phi) - box%cc(1, 1:nc, i_phi)) * &
         box%cc(0, 1:nc, i_eps) / &
         (box%cc(1, 1:nc, i_eps) + box%cc(0, 1:nc, i_eps))
    box%fx(nc+1, :, f_fld) = 2 * inv_dr * &
         (box%cc(nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, i_phi)) * &
         box%cc(nc+1, 1:nc, i_eps) / &
         (box%cc(nc+1, 1:nc, i_eps) + box%cc(nc, 1:nc, i_eps))
    box%fy(:, 1, f_fld) = 2 * inv_dr * &
         (box%cc(1:nc, 0, i_phi) - box%cc(1:nc, 1, i_phi)) * &
         box%cc(1:nc, 0, i_eps) / &
         (box%cc(1:nc, 1, i_eps) + box%cc(1:nc, 0, i_eps))
    box%fy(:, nc+1, f_fld) = 2 * inv_dr * &
         (box%cc(1:nc, nc, i_phi) - box%cc(1:nc, nc+1, i_phi)) * &
         box%cc(1:nc, nc+1, i_eps) / &
         (box%cc(1:nc, nc+1, i_eps) + box%cc(1:nc, nc, i_eps))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fx(1:nc, 1:nc, f_fld) + box%fx(2:nc+1, 1:nc, f_fld))**2 + &
         0.25_dp * (box%fy(1:nc, 1:nc, f_fld) + box%fy(1:nc, 2:nc+1, f_fld))**2)
  end subroutine fld_from_pot

  ! Koren limiter for use in advection scheme
  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, &
         (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  ! Safe division for use in Koren limiter
  elemental function ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: ratio
    if (denominator /= 0.0_dp) then
       ratio = numerator / denominator
    else
       ratio = 1 / epsilon(1.0_dp)
    end if
  end function ratio

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr, theta
    real(dp)                    :: gradp, gradc, gradn
    integer                     :: i, j, nc, nb_id
    type(my_data) :: md

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr
    md     = transfer(boxes(id)%ud, md)

    if (md%box_type == full_diel_box) then
       boxes(id)%fx(:, :, f_elec) = 0
       boxes(id)%fy(:, :, f_elec) = 0
       return
    end if

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1

          gradc = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i-1, j, i_elec)
          if (boxes(id)%fx(i, j, f_fld) * mobility < 0.0_dp) then

             if (i == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hx)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(2, j, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i+1, j, i_elec) - boxes(id)%cc(i, j, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_fld) * mobility * &
                  (boxes(id)%cc(i, j, i_elec) - limiter_koren(theta) * gradn)
          else                  ! boxes(id)%fx(i, j, f_fld) * mobility > 0

             if (i == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_lx)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i-1, j, i_elec) - boxes(nb_id)%cc(nc-1, j, i_elec)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i-1, j, i_elec) - boxes(id)%cc(i-2, j, i_elec)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_fld) * mobility * &
                  (boxes(id)%cc(i-1, j, i_elec) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc

          gradc = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i, j-1, i_elec)

          if (boxes(id)%fy(i, j, f_fld) * mobility < 0.0_dp) then

             if (j == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hy)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(i, 2, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i, j+1, i_elec) - boxes(id)%cc(i, j, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_fld) * mobility * &
                  (boxes(id)%cc(i, j, i_elec) - limiter_koren(theta) * gradn)
          else                  ! boxes(id)%fy(i, j, f_fld) * mobility > 0

             if (j == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_ly)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i, j-1, i_elec) - boxes(nb_id)%cc(i, nc-1, i_elec)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i, j-1, i_elec) - boxes(id)%cc(i, j-2, i_elec)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_fld) * mobility * &
                  (boxes(id)%cc(i, j-1, i_elec) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    if (md%box_type == part_diel_box) then
       do j = 1, nc
          do i = 1, nc
             if (boxes(id)%cc(i, j, i_eps) > 1.0_dp) then
                ! If there's a dielectric, there is no outflux
                if (boxes(id)%fx(i, j, f_elec) < 0) &
                     boxes(id)%fx(i, j, f_elec) = 0
                if (boxes(id)%fx(i+1, j, f_elec) > 0) &
                     boxes(id)%fx(i+1, j, f_elec) = 0
                if (boxes(id)%fy(i, j, f_elec) < 0) &
                     boxes(id)%fy(i, j, f_elec) = 0
                if (boxes(id)%fy(i, j+1, f_elec) > 0) &
                     boxes(id)%fy(i, j+1, f_elec) = 0
             end if
          end do
       end do
    end if

  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, fld
    integer                     :: i, j, nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    do j = 1, nc
       do i = 1, nc
          fld = box%cc(i,j, i_fld)
          src = abs(mobility * fld) * get_alpha(fld) * &
               dt(1) * box%cc(i, j, i_elec)
          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + src
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src
       end do
    end do

    box%cc(1:nc, 1:nc, i_elec) = box%cc(1:nc, 1:nc, i_elec) + dt(1) * ( &
         (box%fx(1:nc, :, f_elec) - box%fx(2:nc+1, :, f_elec)) * inv_dr + &
         (box%fy(:, 1:nc, f_elec) - box%fy(:, 2:nc+1, f_elec)) * inv_dr)

  end subroutine update_solution

  ! Get ionization coefficient (1/m) for a given field strenght
  real(dp) function get_alpha(fld)
    real(dp), intent(in) :: fld
    ! Breakdown fld of 3 MV/m
    get_alpha = max(0.0_dp, 1e5_dp * exp(1 - 1e7_dp/(abs(fld)+epsilon(1.0_dp))) - 9697.2_dp)
  end function get_alpha

  ! For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_children(tree, ref_info)
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, nc

    nc = tree%n_cell

    do lvl = 1, tree%max_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          ! Immediately fill i_eps ghost cells
          call a2_prolong0_to(tree%boxes, id, i_eps, &
               [0,0], [nc+1, nc+1])
          call set_box_type(tree%boxes(id))
          call a2_prolong1_to(tree%boxes, id, i_elec)
          call a2_prolong1_to(tree%boxes, id, i_pion)
          call a2_prolong1_to(tree%boxes, id, i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box_sides(tree%boxes, id, i_elec, &
               a2_sides_interp, sides_bc_dens)
          call a2_gc_box_sides(tree%boxes, id, i_pion, &
               a2_sides_interp, sides_bc_dens)
          call a2_gc_box_sides(tree%boxes, id, i_phi, &
               a2_sides_extrap, sides_bc_pot)
       end do
    end do
  end subroutine prolong_to_new_children

  ! This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       ! Neumann zero
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       ! Neumann zero
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       ! Dirichlet zero
       boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       ! Dirichlet
       boxes(id)%cc(1:nc, nc+1, iv) = 2 * 2.5e6_dp * domain_len &
            - boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc_pot

  ! This fills ghost cells near physical boundaries for the electron density
  subroutine sides_bc_dens(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       ! Neumann zero
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       ! Neumann zero
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       ! Neumann zero
       boxes(id)%cc(1:nc, 0, iv) = boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       ! Neumann zero
       boxes(id)%cc(1:nc, nc+1, iv) = boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc_dens

  ! Based on the box type, apply a Gauss-Seidel relaxation scheme
  subroutine auto_box_gsrb(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: redblack_cntr
    type(mg2_t), intent(in)     :: mg
    type(my_data) :: md

    md = transfer(box%ud, md)
    select case(md%box_type)
    case (full_diel_box, part_diel_box)
       call mg2_box_gsrb_lpld(box, redblack_cntr, mg)
    case (full_gas_box)
       call mg2_box_gsrb_lpl(box, redblack_cntr, mg)
    end select
  end subroutine auto_box_gsrb

  ! Based on the box type, apply the Poisson operator
  subroutine auto_box_op(box, i_out, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_out
    type(mg2_t), intent(in)     :: mg
    type(my_data) :: md

    md = transfer(box%ud, md)
    select case(md%box_type)
    case (full_diel_box, part_diel_box)
       call mg2_box_lpld(box, i_out, mg)
    case (full_gas_box)
       call mg2_box_lpl(box, i_out, mg)
    end select
  end subroutine auto_box_op

  ! Based on the box type, correct the solution of the children
  subroutine auto_box_corr(box_p, box_c, mg)
    type(box2_t), intent(inout) :: box_c
    type(box2_t), intent(in)    :: box_p
    type(mg2_t), intent(in)     :: mg
    type(my_data) :: md

    md = transfer(box_p%ud, md)
    select case(md%box_type)
    case (full_diel_box, part_diel_box)
       call mg2_box_corr_lpld(box_p, box_c, mg)
    case (full_gas_box)
       call mg2_box_corr_lpl(box_p, box_c, mg)
    end select
  end subroutine auto_box_corr

end program test_streamer_2d
