#include "../src/cpp_macros.h"

program KT_euler
  use m_af_all
  implicit none

  integer, parameter :: n_vars = 2+NDIM
  integer, parameter :: ncells = 8
  integer, parameter :: coord_type = af_xyz
  real(dp), parameter :: euler_gamma = 1.4_dp

  integer            :: i_rho, i_mom(NDIM), i_e
  integer            :: variables(n_vars)
  integer            :: fluxes(n_vars)
  real(dp)           :: l_max(NDIM), l_min(NDIM)
  integer            :: grid(NDIM)
  logical            :: periodicBC(NDIM)
  real(dp)           :: u0(4, n_vars)
  type(af_t)         :: tree
  real(dp)           :: dt, time, end_time
  integer            :: t_iter
  character(len=100) :: fname
  integer            :: n
  real(dp)           :: rho_max

  ! AMR stuff
  type(ref_info_t) :: refine_info
  integer          :: refine_steps
  real(dp)         :: dr_min(NDIM)

  print *, "Running Euler 2D with KT scheme"
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "rho", ix=i_rho)
  call af_add_cc_variable(tree, "mom_x", ix=i_mom(1))
  call af_add_cc_variable(tree, "mom_y", ix=i_mom(2))
#if NDIM == 3
  call af_add_cc_variable(tree, "mom_z", ix=i_mom(3))
#endif
  call af_add_cc_variable(tree, "E", ix=i_e)
  variables = [i_rho, i_mom, i_e]

  do n = 1, n_vars
     call af_add_fc_variable(tree, "flux", ix=fluxes(n))
  end do

  call af_set_cc_methods(tree, i_rho, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_mom(1), af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_mom(2), af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_e, af_bc_neumann_zero)

  ! Config 1
  ! u0(:, i_e)      = [1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp]
  ! u0(:, i_rho)    = [1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp]
  ! u0(:, i_mom(1)) = [0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp]
  ! u0(:, i_mom(2)) = [0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp]

  ! Config 6
  ! u0(:, i_e)      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
  ! u0(:, i_rho)    = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
  ! u0(:, i_mom(1)) = [0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp]
  ! u0(:, i_mom(2)) = [-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]

  ! 1D Sod shock test case
  u0(:, i_rho)    = [0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp]
  u0(:, i_e)      = [0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp]
  u0(:, i_mom(1)) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  u0(:, i_mom(2)) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

  call to_conservative(4, n_vars, u0)

  grid(:)       = 50*ncells
  l_max(:)      = 1.0_dp
  l_min(:)      = 0.0_dp
  periodicBC(:) = .false.

  call af_init(tree, ncells, l_max-l_min, grid, &
       periodic=periodicBC, r_min=l_min, &
       coord=coord_type)

  !Init mesh refinement
  !  do
  !      refine_steps = refine_steps + 1
  !      !Settng init conds for each refinement is needed as we use that data as
  !      !refinement criterion
  !      call af_loop_box(tree, set_init_conds)
  !      call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])
  !      call af_adjust_refinement(tree, ref_rout, refine_info, 1)
  !  if (refine_info%n_add == 0) exit
  !  end do
  !  call af_restrict_tree(tree, i_rho)
  !  call af_restrict_tree(tree, i_mom(1))
  !  call af_restrict_tree(tree, i_mom(2))
  !  call af_restrict_tree(tree, i_e)
  !  call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])


  call af_loop_box(tree, set_init_conds)

  call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])

  !call af_write_silo(tree, 'eulerInit', dir='output')

  !Setting the timestep data
  time = 0.0_dp
  end_time = 0.2_dp
  t_iter = 0
  dt = 2.0e-04_dp
  do
     if (mod(t_iter, 10) == 0) then
        write(fname, "(A,I0)") "KT_euler_" // DIMNAME // "_", t_iter
        call af_write_silo(tree, trim(fname), t_iter, time, dir="output", &
             add_vars = write_primitives, add_names=["xVel","yVel","pres"])
     end if

     call af_loop_tree(tree, compute_flux)
     call af_consistent_fluxes(tree, variables)
     call af_loop_box_arg(tree, update_solution, [dt])
     do n = 1, n_vars
        call af_restrict_tree(tree, variables(n))
     end do
     call af_gc_tree(tree, variables)
     !call af_adjust_refinement(tree, ref_rout, refine_info, 1)

     !call af_loop_box_arg(tree, update_solution, [dt])
     !call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])

     t_iter = t_iter + 1
     time = time + dt

     call af_tree_maxabs_cc(tree, i_rho, rho_max)
     if (rho_max > 10.0_dp) &
          error stop "solution diverging!"

     if (time > end_time) exit
  end do

  call af_destroy(tree)

contains

  subroutine set_init_conds(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0, nc+1)
       rr = af_r_cc(box, [IJK])
       if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
          box%cc(IJK, i_rho)    = u0(1, 1)
          box%cc(IJK, i_mom(1)) = u0(1, 2)
          box%cc(IJK, i_mom(2)) = u0(1, 3)
          box%cc(IJK, i_e)      = u0(1, 4)
       elseif (rr(1) <= 0.5_dp .and. rr(2) >= 0.5_dp) then
          box%cc(IJK, i_rho)    = u0(2, 1)
          box%cc(IJK, i_mom(1)) = u0(2, 2)
          box%cc(IJK, i_mom(2)) = u0(2, 3)
          box%cc(IJK, i_e)      = u0(2, 4)
       elseif (rr(1) <= 0.5_dp .and. rr(2) <= 0.5_dp) then
          box%cc(IJK, i_rho)    = u0(3, 1)
          box%cc(IJK, i_mom(1)) = u0(3, 2)
          box%cc(IJK, i_mom(2)) = u0(3, 3)
          box%cc(IJK, i_e)      = u0(3, 4)
       else
          box%cc(IJK, i_rho)    = u0(4, 1)
          box%cc(IJK, i_mom(1)) = u0(4, 2)
          box%cc(IJK, i_mom(2)) = u0(4, 3)
          box%cc(IJK, i_e)      = u0(4, 4)
       end if
    end do; CLOSE_DO
  end subroutine set_init_conds

  subroutine compute_flux(tree, id)
    use m_af_flux_schemes
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id

    call flux_generic(tree, id, tree%n_cell, n_vars, variables, fluxes, &
         to_primitive, to_conservative, get_max_wavespeed_1d, get_fluxes_1d)
  end subroutine compute_flux

  subroutine to_primitive(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)

    u(:, i_mom(1)) = u(:, i_mom(1))/u(:, i_rho)
    u(:, i_mom(2)) = u(:, i_mom(2))/u(:, i_rho)
    u(:, i_e) = (euler_gamma-1.0_dp) * (u(:, i_e) - &
         0.5_dp*u(:, i_rho)* sum(u(:, i_mom(:))**2, dim=2))
  end subroutine to_primitive

  subroutine to_conservative(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
    real(dp)                :: kin_en(n_values)
    real(dp)                :: inv_fac
    integer                 :: i

    ! Compute kinetic energy (0.5 * rho * velocity^2)
    kin_en = 0.5_dp * u(:, i_rho) * sum(u(:, i_mom(:))**2, dim=2)

    ! Compute energy from pressure and kinetic energy
    inv_fac = 1/(euler_gamma - 1.0_dp)
    u(:, i_e) = u(:, i_e) * inv_fac + kin_en

    ! Compute momentum from density and velocity components
    do i = 1, NDIM
       u(:, i_mom(i)) = u(:, i_rho) * u(:, i_mom(i))
    end do
  end subroutine to_conservative

  subroutine get_max_wavespeed_1d(n_values, n_var, flux_dim, u, w)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var) !< Primitive variables
    real(dp), intent(out) :: w(n_values) !< Maximum speed
    real(dp)              :: sound_speeds(n_values)

    sound_speeds = sqrt(euler_gamma * u(:, i_e) / u(:, i_rho))
    w = sound_speeds + abs(u(:, i_mom(flux_dim)))
  end subroutine get_max_wavespeed_1d

  subroutine get_fluxes_1d(n_values, n_var, flux_dim, u, flux)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var)
    real(dp), intent(out) :: flux(n_values, n_var)
    real(dp)              :: E(n_values), inv_fac
    integer               :: i

    ! Compute left and right flux for conservative variables from the primitive
    ! reconstructed values.

    ! Density flux
    flux(:, i_rho) = u(:, i_rho) * u(:, i_mom(flux_dim))

    ! Momentum flux
    do i = 1, NDIM
       flux(:,  i_mom(i)) = u(:, i_rho) * &
            u(:, i_mom(i)) * u(:, i_mom(flux_dim))
    end do

    ! Add pressure term
    flux(:, i_mom(flux_dim)) = flux(:, i_mom(flux_dim)) + u(:, i_e)

    ! Compute energy
    inv_fac = 1/(euler_gamma-1.0_dp)
    E = u(:, i_e) * inv_fac + 0.5_dp * u(:, i_rho) * &
         sum(u(:, i_mom(:))**2, dim=2)

    ! Energy flux
    flux(:, i_e) = u(:, i_mom(flux_dim)) * (E + u(:, i_e))

  end subroutine get_fluxes_1d

  subroutine update_solution(box, dt)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt(:)
    real(dp)                   :: inv_dr(NDIM)
    integer                    :: IJK, nc, n, iv, iflux

    nc     = box%n_cell
    inv_dr = 1.0_dp/box%dr

    do n = 1, n_vars
       iv    = variables(n)
       iflux = fluxes(n)

       do KJI_DO(1, nc)
#if NDIM == 2
          box%cc(IJK, iv) = box%cc(IJK, iv) - dt(1) * ( &
               inv_dr(1) * &
               (box%fc(i+1, j, 1, iflux) - box%fc(i, j, 1, iflux)) + &
               inv_dr(2) * &
               (box%fc(i, j+1, 2, iflux) - box%fc(i, j, 2, iflux)))
#elif NDIM == 3
          box%cc(IJK, iv) = box%cc(IJK, iv) - dt(1) * ( &
               inv_dr(1) * &
               (box%fc(i+1, j, k, 1, iflux) - box%fc(i, j, k, 1, iflux)) + &
               inv_dr(2) * &
               (box%fc(i, j+1, k, 2, iflux) - box%fc(i, j, k, 2, iflux)) + &
               inv_dr(3) * &
               (box%fc(i, j, k+1, 3, iflux) - box%fc(i, j, k, 3, iflux)))
#endif
       end do; CLOSE_DO
    end do
  end subroutine update_solution

  subroutine ref_rout( box, cell_flags )
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: diff, tol
    integer                 :: IJK, nc

    nc  = box%n_cell
    tol = 1.0e-6_dp
    do KJI_DO(1,nc)
       diff =  box%dr(1)**2*abs(box%cc(i+1,j,i_rho)+box%cc(i-1,j,i_rho) &
            -2*box%cc(i,j,i_rho)) + &
            box%dr(2)**2*abs(box%cc(i,j+1,i_rho)+box%cc(i,j-1,i_rho) &
            -2*box%cc(i,j,i_rho))
       if (diff > tol .and. box%lvl .le. 4) then
          cell_flags(IJK) = af_do_ref
       else if (diff < 0.1_dp*tol) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO
  end subroutine ref_rout

  subroutine write_primitives(box, new_vars, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: n_var
    real(dp)                :: new_vars(DTIMES(0:box%n_cell+1),n_var)

    ! X Velocity
    new_vars(DTIMES(:), 1) = box%cc(DTIMES(:), i_mom(1))/box%cc(DTIMES(:), i_rho)
    ! Y Velocity
    new_vars(DTIMES(:), 2) = box%cc(DTIMES(:), i_mom(2))/box%cc(DTIMES(:), i_rho)
    ! Pressure
    new_vars(DTIMES(:), 3) = (euler_gamma-1.0_dp)*(box%cc(DTIMES(:), i_e) - &
         sum(box%cc(DTIMES(:), i_mom(:))**2, dim=NDIM+1) &
         /(2.0_dp*box%cc(DTIMES(:), i_rho)))

  end subroutine write_primitives


end program KT_euler
