!> \example simple_streamer_2d.f90
!>
!> A simplified model for ionization waves and/or streamers in 2D
program simple_streamer_2d

  use m_af_types
  use m_af_core
  use m_af_ghostcell
  use m_af_utils
  use m_af_restrict
  use m_af_multigrid
  use m_af_output
  use m_write_silo
  use m_af_prolong

  implicit none

  integer            :: i, n
  character(len=100) :: fname
  type(af_t)         :: tree ! This contains the full grid information
  type(mg_t)         :: mg   ! Multigrid option struct
  type(ref_info_t)   :: refine_info

  ! Indices of cell-centered variables
  integer :: i_elec     ! Electron density
  integer :: i_pion     ! Positive ion density
  integer :: i_elec_old ! For time-stepping scheme
  integer :: i_pion_old ! For time-stepping scheme
  integer :: i_phi      ! Electrical potential
  integer :: i_fld      ! Electric field norm
  integer :: i_rhs      ! Source term Poisson

  ! Indices of face-centered variables **
  integer :: f_elec ! Electron flux
  integer :: f_fld  ! Electric field vector

  ! Simulation parameters
  real(dp), parameter :: end_time      = 3e-9_dp
  real(dp), parameter :: dt_output     = 20e-11_dp
  real(dp), parameter :: dt_max        = 20e-11_dp
  integer, parameter  :: ref_per_steps = 2
  integer, parameter  :: box_size      = 8

  ! Physical parameters
  real(dp), parameter :: applied_field = -0.8e7_dp
  real(dp), parameter :: mobility      = 0.03_dp
  real(dp), parameter :: diffusion_c   = 0.2_dp

  ! Computational domain
  real(dp), parameter :: domain_length = 2e-3_dp
  real(dp), parameter :: refine_max_dx = 1e-3_dp
  real(dp), parameter :: refine_min_dx = 1e-9_dp

  ! Settings for the initial conditions
  real(dp), parameter :: init_density = 1e15_dp
  real(dp), parameter :: init_y_min   = 0.8_dp * domain_length
  real(dp), parameter :: init_y_max   = 0.9_dp * domain_length

  ! Simulation variables
  real(dp) :: dt
  real(dp) :: time
  integer  :: output_count

  ! To test charge conservation
  real(dp) :: sum_elec, sum_pion

  call af_add_cc_variable(tree, "elec", ix=i_elec, n_copies=2)
  i_elec_old = af_find_cc_variable(tree, "elec_2")
  call af_add_cc_variable(tree, "pion", ix=i_pion, n_copies=2)
  i_pion_old = af_find_cc_variable(tree, "pion_2")
  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "fld", ix=i_fld)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)

  call af_add_fc_variable(tree, "f_elec", ix=f_elec)
  call af_add_fc_variable(tree, "f_fld", ix=f_fld)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs

  ! Routines to use for ...
  mg%sides_bc => sides_bc_pot ! Filling ghost cell on physical boundaries

  ! This routine always needs to be called when using multigrid
  call mg_init_mg(mg)

  call af_set_cc_methods(tree, i_elec, af_bc_dirichlet_zero, &
       af_gc_interp_lim, af_prolong_limit)
  call af_set_cc_methods(tree, i_pion, af_bc_dirichlet_zero, &
       af_gc_interp_lim, af_prolong_limit)
  call af_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)

  output_count = 0 ! Number of output files written
  time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     ! For each box in tree, set the initial conditions
     call af_loop_box(tree, set_initial_condition)

     ! Compute electric field on the tree.
     ! First perform multigrid to get electric potential,
     ! then take numerical gradient to geld field.
     call compute_fld(tree, .false.)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine af_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     call af_adjust_refinement(tree, &               ! tree
          refinement_routine, & ! Refinement function
          refine_info)          ! Information about refinement

     ! If no new boxes have been added or removed, exit the loop
     if (refine_info%n_add == 0 .and. refine_info%n_rm == 0) exit
  end do

  call af_print_info(tree)

  do
     ! Get a new time step, which is at most dt_max
     call af_reduction(tree, &    ! Tree to do the reduction on
          max_dt, &  ! function
          get_min, & ! function
          dt_max, &  ! Initial value for the reduction
          dt)        ! Result of the reduction

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?", dt
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_count * dt_output <= time) then
        output_count = output_count + 1
        write(fname, "(A,I6.6)") "simple_streamer_2d_", output_count

        ! Write the cell centered data of a tree to a Silo file. Only the
        ! leaves of the tree are used
        call af_write_silo(tree, fname, output_count, time, dir="output")

        call af_tree_sum_cc(tree, i_elec, sum_elec)
        call af_tree_sum_cc(tree, i_pion, sum_pion)
        print *, "sum(pion-elec)/sum(pion)", (sum_pion - sum_elec)/sum_pion
     end if

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, ref_per_steps
        time = time + dt

        ! Copy previous solution
        call af_tree_copy_cc(tree, i_elec, i_elec_old)
        call af_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! First calculate fluxes
           call af_loop_boxes(tree, fluxes_koren, leaves_only=.true.)
           call af_consistent_fluxes(tree, [f_elec])

           ! Update the solution
           call af_loop_box_arg(tree, update_solution, [dt], &
                leaves_only=.true.)

           ! Compute new field on first iteration
           if (i == 1) call compute_fld(tree, .true.)
        end do

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call af_loop_box(tree, average_dens)

        ! Compute field with new density
        call compute_fld(tree, .true.)
     end do

     ! Fill ghost cells for i_pion and i_elec
     call af_gc_tree(tree, i_elec)
     call af_gc_tree(tree, i_pion)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine af_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call af_adjust_refinement(tree, &               ! tree
          refinement_routine, & ! Refinement function
          refine_info, &        ! Information about refinement
          4)                    ! Buffer width (in cells)

     if (refine_info%n_add > 0 .or. refine_info%n_rm > 0) then
        ! Compute the field on the new mesh
        call compute_fld(tree, .true.)
     end if
  end do

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(af_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = domain_length / box_size

    ! Initialize tree
    call af_init(tree, &                   ! The tree to initialize
         box_size, &               ! Boxes have box_size^dim cells
         dr, &                     ! spacing of a cell at lvl 1
         coarsen_to=2, &           ! Create additional coarse grids down to this size.
         n_boxes = n_boxes_init)   ! Allocate initial storage for n_boxes.

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...

    nb_list(af_neighb_lowy, id)  = -1 ! physical boundary
    nb_list(af_neighb_highy, id) = -1 ! idem
    nb_list(af_neighb_lowx, id)  = id ! periodic boundary
    nb_list(af_neighb_highx, id) = id ! idem

    ! Create the base mesh
    call af_set_base(tree, 1, ix_list, nb_list)

  end subroutine init_tree

  !> This routine sets the refinement flag for boxes(id)
  subroutine refinement_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: dx, dens, fld, adx, xy(2)

    nc      = box%n_cell
    dx      = box%dr

    do j = 1, nc
       do i = 1, nc
          xy   = af_r_cc(box, [i,j])
          dens = box%cc(i, j, i_elec)
          fld = box%cc(i, j, i_fld)
          adx = get_alpha(fld) * dx

          if (dens > 1e0_dp .and. adx > 0.8_dp) then
             cell_flags(i, j) = af_do_ref
          else if (dx < 1.25e-5_dp .and. adx < 0.1_dp) then
             cell_flags(i, j) = af_rm_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    else if (dx > 0.5_dp * refine_max_dx) then
       where(cell_flags == af_rm_ref) cell_flags = af_keep_ref
    end if

  end subroutine refinement_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), normal_rands(2)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = af_r_cc(box, [i,j])

          if (xy(2) > init_y_min .and. xy(2) < init_y_max) then
             ! Approximate Poisson distribution with normal distribution
             normal_rands = two_normals(box%dr**3 * init_density, &
                  sqrt(box%dr**3 * init_density))
             ! Prevent negative numbers
             box%cc(i, j, i_elec) = abs(normal_rands(1)) / box%dr**3
          else
             box%cc(i, j, i_elec) = 0
          end if
       end do
    end do

    box%cc(:, :, i_pion) = box%cc(:, :, i_elec)
    box%cc(:, :, i_phi)  = 0 ! Inital potential set to zero

  end subroutine set_initial_condition

  !> Return two normal random variates
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function two_normals(mean, sigma) result(rands)
    real(dp), intent(in) :: mean, sigma
    real(dp) :: rands(2), sum_sq

    do
       call random_number(rands)
       rands = 2 * rands - 1
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
    rands = mean + rands * sigma
  end function two_normals

  !> This function computes the minimum val of a and b
  real(dp) function get_min(a, b)
    real(dp), intent(in) :: a, b
    get_min = min(a, b)
  end function get_min

  !> Get maximum time step based on e.g. CFL criteria
  real(dp) function max_dt(box)
    type(box_t), intent(in) :: box
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    integer :: i, j, nc
    real(dp)               :: fld(2)
    real(dp)               :: dt_cfl, dt_dif, dt_drt

    nc = box%n_cell
    dt_cfl = dt_max

    do j = 1, nc
       do i = 1, nc
          fld(1) = 0.5_dp * (box%fc(i, j, 1, f_fld) + &
               box%fc(i+1, j, 1, f_fld))
          fld(2) = 0.5_dp * (box%fc(i, j, 2, f_fld) + &
               box%fc(i, j+1, 2, f_fld))

          ! The 0.5 is here because of the explicit trapezoidal rule
          dt_cfl = min(dt_cfl, 0.5_dp / sum(abs(fld * mobility) / box%dr))
       end do
    end do

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * mobility * &
         maxval(box%cc(1:nc, 1:nc, i_elec)) + epsilon(1.0_dp))

    ! Diffusion condition
    dt_dif = 0.25_dp * box%dr**2 / max(diffusion_c, epsilon(1.0_dp))

    max_dt = min(dt_cfl, dt_drt, dt_dif)
  end function max_dt


  !> This function gets the alpha value
  !>
  !> Taken from: Spatially hybrid computations for streamer discharges: II. Fully
  !> 3D simulations, Chao Li, Ute Ebert, Willem Hundsdorfer, J. Comput. Phys.
  !> 231, 1020-1050 (2012), doi:10.1016/j.jcp.2011.07.023
  elemental function get_alpha(fld) result(alpha)
    real(dp), intent(in) :: fld
    real(dp)             :: alpha
    real(dp), parameter  :: c0 = 1.04e1_dp
    real(dp), parameter  :: c1 = 0.601_dp
    real(dp), parameter  :: c2 = 1.86e7_dp

    alpha = exp(c0) * (abs(fld) * 1e-5_dp)**c1 * exp(-c2 / abs(fld))
  end function get_alpha

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, have_guess)
    type(af_t), intent(inout) :: tree
    logical, intent(in)       :: have_guess
    real(dp), parameter       :: fac = 1.6021766208e-19_dp / 8.8541878176e-12_dp
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
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

    ! Perform an FMG cycle (logicals: store residual, first call)
    call mg_fas_fmg(tree, mg, .false., have_guess)

    ! Compute field from potential
    call af_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call af_gc_tree(tree, i_fld, af_gc_interp, af_bc_neumann_zero)
  end subroutine compute_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fc(1:nc+1, 1:nc, 1, f_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, f_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fc(1:nc, 1:nc, 1, f_fld) + &
         box%fc(2:nc+1, 1:nc, 1, f_fld))**2 + &
         0.25_dp * (box%fc(1:nc, 1:nc, 2, f_fld) + &
         box%fc(1:nc, 2:nc+1, 2, f_fld))**2)
  end subroutine fld_from_pot

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    use m_flux_schemes
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr
    real(dp)                    :: cc(-1:boxes(id)%n_cell+2, -1:boxes(id)%n_cell+2)
    real(dp), allocatable       :: v(:, :, :), dc(:, :, :)
    integer                     :: nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    allocate(v(1:nc+1, 1:nc+1, 2))
    allocate(dc(1:nc+1, 1:nc+1, 2))

    call af_gc_box(boxes, id, i_elec, af_gc_interp_lim, &
         af_bc_dirichlet_zero)
    call af_gc2_box(boxes, id, i_elec, af_gc2_prolong_linear, &
         af_bc2_dirichlet_zero, cc, nc)

    v = -mobility * boxes(id)%fc(:, :, :, f_fld)
    dc = diffusion_c

    call flux_koren_2d(cc, v, nc, 2)
    call flux_diff_2d(cc, dc, inv_dr, nc, 2)

    boxes(id)%fc(:, :, :, f_elec) = v + dc
  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, sflux, fld
    real(dp)                    :: alpha
    integer                     :: i, j, nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    do j = 1, nc
       do i = 1, nc
          fld   = box%cc(i,j, i_fld)
          alpha = get_alpha(fld)
          src   = box%cc(i, j, i_elec) * mobility * abs(fld) * alpha

          sflux = (sum(box%fc(i, j, :, f_elec)) - &
               box%fc(i+1, j, 1, f_elec) - &
               box%fc(i, j+1, 2, f_elec)) * inv_dr

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + (src + sflux) * dt(1)
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src * dt(1)
       end do
    end do

  end subroutine update_solution

  !> This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(box, nb, iv, bc_type)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (af_neighb_lowy)
       bc_type = af_bc_neumann
       box%cc(1:nc, 0, iv) = applied_field
    case (af_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = 0
    case default
       stop "sides_bc_pot: unspecified boundary"
    end select
  end subroutine sides_bc_pot

end program simple_streamer_2d
