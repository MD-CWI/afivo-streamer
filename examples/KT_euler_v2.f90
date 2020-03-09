#include "../src/cpp_macros.h"

program KT_euler
  use m_af_all
  implicit none

  integer, parameter :: n_vars   = 4
  integer, parameter :: n_gc = 2
  !integer, parameter :: i_rho    = 1
  !integer, parameter :: i_mom(2) = [2,3]
  !integer, parameter :: i_e      = 4
  integer :: if1, if2, if3, if4
  integer :: i_rho, i_mom(2), i_e

  integer, parameter :: ncells = 8
  integer, parameter :: coord_type = af_xyz
  real(dp) :: l_max(NDIM), l_min(NDIM)
  integer :: grid(NDIM)
  logical :: periodicBC(NDIM)
  real(dp), parameter :: Y = 1.4_dp
  real(dp):: p(4), rho(4), u(4), v(4)
  real(dp) :: wSp, error


  type(af_t) :: tree
  real(dp) :: dt, time, end_time
  integer :: t_iter
  character(len=100) :: fname

  !AMR stuff
  type(ref_info_t) :: refine_info
  integer :: refine_steps
  real(dp) :: dr_min(NDIM)

  print *, "Running Euler 2D with KT scheme"
  print *, "Number of threads", af_get_max_threads()

  !Config 1
!  p = (/1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp/)
!  rho = (/1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp/)
!  u = (/0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp/)
!  v = (/0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp/)

  ! Config 6
!  p = (/1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/)
!  rho = (/1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp/)
!  u = (/0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp/)
!  v = (/-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp/)

   !1D Sod shock test case
  rho = (/0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp/)
  p = (/0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp/)
  u = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
  v = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)


  wSp = calc_speed(rho, u, v, p)



  grid(:) = 50*ncells
  l_max(:) = 1.0_dp
  l_min(:) = 0.0_dp
  periodicBC(:) = .false.

  call af_add_cc_variable(tree, "rho", ix=i_rho)
  call af_add_cc_variable(tree, "rhoU", ix=i_mom(1))
  call af_add_cc_variable(tree, "rhoV", ix=i_mom(2))
  call af_add_cc_variable(tree, "E", ix=i_e)

  call af_add_fc_variable(tree, "flux1", ix=if1)
  call af_add_fc_variable(tree, "flux2", ix=if2)
  call af_add_fc_variable(tree, "flux3", ix=if3)
  call af_add_fc_variable(tree, "flux4", ix=if4)

  call af_set_cc_methods(tree, i_rho, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_mom(1), af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_mom(2), af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_e, af_bc_neumann_zero)

  call af_init(tree, ncells, l_max, grid, &
               periodic = periodicBC, r_min= l_min, &
               coord = coord_type)


  !Init mesh refinement
!  do
!      refine_steps = refine_steps + 1
!      !Settng init conds for each refinement is needed as we use that data as
!      !refinement criterion
!      call af_loop_box(tree, setInitConds)
!      call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])
!      call af_adjust_refinement(tree, ref_rout, refine_info, 1)
!  if (refine_info%n_add == 0) exit
!  end do
!  call af_restrict_tree(tree, i_rho)
!  call af_restrict_tree(tree, i_mom(1))
!  call af_restrict_tree(tree, i_mom(2))
!  call af_restrict_tree(tree, i_e)
!  call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])


  call af_loop_box(tree, setInitConds)

  call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])

  !call af_write_silo(tree, 'eulerInit', dir='output')

  !Setting the timestep data
  time = 0.0_dp
  end_time = 0.2_dp
  t_iter = 0
  dt = 2.0e-04_dp
  do
   !Updating the primitive vars
   !call af_tree_copy_cc(tree, i_mom(1), ip2)
   !call af_tree_copy_cc(tree, i_mom(2), ip3)
   !call af_tree_copy_cc(tree, i_e, ip4)
   !call af_loop_box(tree, updatePrimitives)
   !dr_min = af_min_dr(tree)
   !dt = end_time/(sum(wSp/dr_min) + epsilon(1.0_dp))
   !print *, dt
    if (mod(t_iter, 10) == 0) then
      write(fname, "(A,I0)") "KT_euler_" // DIMNAME // "_", t_iter
      call af_write_silo(tree, trim(fname), t_iter, time, dir="output", &
           add_vars = writePrimitives, add_names=["xVel","yVel","pres"])
    end if

   call af_loop_tree(tree, fluxComputation)
   call af_consistent_fluxes(tree, [i_rho,i_mom(1), i_mom(2), i_e])
   call af_loop_box_arg(tree, updateSoln, [dt])
   call af_restrict_tree(tree, i_rho)
   call af_restrict_tree(tree, i_mom(1))
   call af_restrict_tree(tree, i_mom(2))
   call af_restrict_tree(tree, i_e)
   call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])
   !call af_adjust_refinement(tree, ref_rout, refine_info, 1)

   !call af_loop_box_arg(tree, updateSoln, [dt])
   !call af_gc_tree(tree, [i_rho, i_mom(1), i_mom(2), i_e])

   t_iter = t_iter + 1
   time = time + dt


   call af_tree_maxabs_cc(tree, i_rho, error)
   !print *, "Max value of fc: ", test
   if (error > 10.0_dp) then
    print *, "solution diverging!"
    exit
   end if

   if (time > end_time) exit
  end do




  call af_destroy(tree)
  !=====================================================================
  contains


  function calc_speed(rho, u, v, p) result(waveSpeed)
    real(dp), intent(in) :: rho(4), u(4), v(4), p(4)
    real(dp):: waveSpeed
    integer :: i

    waveSpeed = 0.0_dp
    do i=1,4
      waveSpeed = max(waveSpeed, abs(sqrt((Y*p(i))/(rho(i))) + &
                                        sqrt(u(i)**2 + v(i)**2)))
    end do

  end function calc_speed
  !=====================================================================
  function convert_to_conservatives( rho, u, v, p ) result(c)
    real(dp), intent(in) :: rho, u, v, p
    real(dp):: c(4)

    c(1) = rho
    c(2) = rho*u
    c(3) = rho*v
    c(4) = (p/(Y-1.0_dp)) + 0.5_dp*rho*(u**2 + v**2)

  end function convert_to_conservatives

  subroutine setInitConds( box )
    type(box_t), intent(inout) :: box
    integer :: IJK, nc
    real(dp) :: rr(NDIM)
    real(dp) :: conservatives(4)


    nc = box%n_cell
    do KJI_DO(0, nc+1)
      rr = af_r_cc(box, [IJK])
      if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
        conservatives = convert_to_conservatives(rho(1), u(1), v(1), p(1))
        box%cc(IJK, i_rho) = conservatives(1)
        box%cc(IJK, i_mom(1)) = conservatives(2)
        box%cc(IJK, i_mom(2)) = conservatives(3)
        box%cc(IJK, i_e) = conservatives(4)
!        box%cc(IJK, ip2) = u(1)
!        box%cc(IJK, ip3) = v(1)
!        box%cc(IJK, ip4) = p(1)
      elseif (rr(1) .le. 0.5_dp .and. rr(2) .ge. 0.5_dp) then
        conservatives = convert_to_conservatives(rho(2), u(2), v(2), p(2))
        box%cc(IJK, i_rho) = conservatives(1)
        box%cc(IJK, i_mom(1)) = conservatives(2)
        box%cc(IJK, i_mom(2)) = conservatives(3)
        box%cc(IJK, i_e) = conservatives(4)
!        box%cc(IJK, ip2) = u(2)
!        box%cc(IJK, ip3) = v(2)
!        box%cc(IJK, ip4) = p(2)
      elseif (rr(1) .le. 0.5_dp .and. rr(2) .le. 0.5_dp) then
        conservatives = convert_to_conservatives(rho(3), u(3), v(3), p(3))
        box%cc(IJK, i_rho) = conservatives(1)
        box%cc(IJK, i_mom(1)) = conservatives(2)
        box%cc(IJK, i_mom(2)) = conservatives(3)
        box%cc(IJK, i_e) = conservatives(4)
!        box%cc(IJK, ip2) = u(3)
!        box%cc(IJK, ip3) = v(3)
!        box%cc(IJK, ip4) = p(3)
      else
        conservatives = convert_to_conservatives(rho(4), u(4), v(4), p(4))
        box%cc(IJK, i_rho) = conservatives(1)
        box%cc(IJK, i_mom(1)) = conservatives(2)
        box%cc(IJK, i_mom(2)) = conservatives(3)
        box%cc(IJK, i_e) = conservatives(4)
!        box%cc(IJK, ip2) = u(4)
!        box%cc(IJK, ip3) = v(4)
!        box%cc(IJK, ip4) = p(4)
      end if
    end do; CLOSE_DO
  end subroutine setInitConds
  !=====================================================================
  subroutine fluxComputation( tree, id )
    use m_af_flux_schemes
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer                   :: nc, i, j
    real(dp), allocatable     :: cc(DTIMES(:), :)
    real(dp), allocatable     :: prim_vars(:,:), cc_line(:, :)
    real(dp), allocatable     :: u_lr(:,:,:), w_lr(:)
    real(dp), allocatable     :: flux_lr(:,:,:), flux(:,:)


    nc = tree%boxes(id)%n_cell
    allocate(cc(DTIMES(-1:nc+2), n_vars))
    allocate(cc_line(-1:nc+2, n_vars))
    allocate(prim_vars(-1:nc+2, n_vars))
    allocate(u_lr(1:nc+1, 2, n_vars))
    allocate(w_lr(1:nc+1))
    allocate(flux_lr(1:nc+1, 2, n_vars))
    allocate(flux(1:nc+1, n_vars))

    call af_gc2_box(tree, id, [i_rho, i_mom(1), i_mom(2), i_e], cc)

    ! x-direction
    do j = 1, nc
       !We need primitives at the ghost cells as they're used for reconstruction
       cc_line = cc(:, j, :)
       call to_primitive(nc+2*n_gc, n_vars, cc_line, prim_vars)
       call reconstruct_lr_1d(nc, n_gc, n_vars, prim_vars, u_lr)
       call get_max_wavespeed_lr_1d(nc+1, n_vars, 1, u_lr, w_lr)
       call get_fluxes_lr_1d(nc+1, n_vars, 1, u_lr, flux_lr)
       call to_conservatives(nc+1, n_vars, u_lr)
       call flux_kurganovTadmor_1d(nc+1, n_vars, flux_lr, u_lr, w_lr, flux)
       tree%boxes(id)%fc(:, j, 1, :) = flux
    end do

    ! y-direction
    do i = 1, nc
       cc_line = cc(i, :, :)
       call to_primitive(nc+2*n_gc, n_vars, cc_line, prim_vars)
       call reconstruct_lr_1d(nc, n_gc, n_vars, prim_vars, u_lr)
       call get_max_wavespeed_lr_1d(nc+1, n_vars, 2, u_lr, w_lr)
       call get_fluxes_lr_1d(nc+1, n_vars, 2, u_lr, flux_lr)
       call to_conservatives(nc+1, n_vars, u_lr)
       call flux_kurganovTadmor_1d(nc+1, n_vars, flux_lr, u_lr, w_lr, flux)
       tree%boxes(id)%fc(i, :, 2, :) = flux
    end do

  end subroutine fluxComputation

  subroutine to_primitive(n_values, n_vars, cc, prim_vars)
    integer, intent(in) :: n_values, n_vars
    real(dp), intent(in) :: cc(n_values, n_vars)
    real(dp), intent(out) :: prim_vars(n_values, n_vars)

    prim_vars(:, i_rho) = cc(:, i_rho)
    prim_vars(:, i_mom(1)) = cc(:, i_mom(1))/cc(:, i_rho)
    prim_vars(:, i_mom(2)) = cc(:, i_mom(2))/cc(:, i_rho)
    prim_vars(:, i_e) = (Y-1.0_dp) * (cc(:, i_e) - &
         0.5_dp*prim_vars(:, i_rho)* sum(prim_vars(:, i_mom(:))**2, dim=2))
  end subroutine to_primitive

  subroutine to_conservatives(nf, n_vars, u_lr)
    integer, intent(in)     :: nf, n_vars
    real(dp), intent(inout) :: u_lr(nf, 2, n_vars)
    real(dp)                :: kin_en(nf, 2)
    integer                 :: i

    ! Compute 0.5 rho velocity^2
    kin_en = 0.5_dp * u_lr(:, :, i_rho) * sum(u_lr(:, :, i_mom(:))**2, dim=3)

    ! Compute energy from pressure and kinetic energy
    u_lr(:,:, i_e) = u_lr(:,:, i_e)/(Y - 1.0_dp) + kin_en

    ! Compute momentum from density and velocity components
    do i = 1, NDIM
       u_lr(:,:, i_mom(i)) = u_lr(:,:,i_rho)*u_lr(:,:, i_mom(i))
    end do
  end subroutine to_conservatives

  subroutine get_max_wavespeed_lr_1d(nf, n_var, flux_dim, u_lr, w_lr)
    integer, intent(in)   :: nf    !< Number of cell faces
    integer, intent(in)   :: n_var !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u_lr(nf, 2, n_var)
    real(dp), intent(out) :: w_lr(nf)
    integer               :: n
    real(dp)              :: sound_speeds(nf, 2)

    sound_speeds = sqrt(Y * u_lr(:, :, i_e) / u_lr(:, :, i_rho))
    w_lr = maxval(sound_speeds + abs(u_lr(:, :, i_mom(flux_dim))),2)
  end subroutine get_max_wavespeed_lr_1d

  subroutine get_fluxes_lr_1d(nf, n_var, flux_dim, u_lr, flux_lr)
    integer, intent(in)   :: nf       !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u_lr(nf, 2, n_var)
    real(dp), intent(out) :: flux_lr(nf, 2, n_var)
    real(dp)              :: E(nf, 2)
    integer               :: i

    ! Compute left and right flux for conservative variables from the primitive
    ! reconstructed values.

    ! Density flux
    flux_lr(:,:,i_rho) = u_lr(:,:,i_rho)*u_lr(:,:,i_mom(flux_dim))

    ! Momentum flux
    do i = 1, NDIM
       flux_lr(:,:, i_mom(i)) = u_lr(:,:,i_rho) * &
            u_lr(:,:,i_mom(i)) * u_lr(:,:,i_mom(flux_dim))
    end do

    ! Add pressure term
    flux_lr(:,:, i_mom(flux_dim)) = flux_lr(:,:, i_mom(flux_dim)) + u_lr(:,:,i_e)

    ! Compute energy
    E = u_lr(:,:,i_e)/(Y-1.0_dp) + 0.5_dp*u_lr(:,:,i_rho)*(&
         u_lr(:,:,i_mom(1))**2 + u_lr(:,:,i_mom(2))**2)

    ! Energy flux
    flux_lr(:,:, i_e) = u_lr(:,:,i_mom(flux_dim))*(E + u_lr(:,:,i_e))

  end subroutine get_fluxes_lr_1d

  subroutine updateSoln( box, dt )
    type(box_t), intent(inout) :: box
    real(dp), intent(in) :: dt(:)
    real(dp) :: inv_dr(NDIM), avg
    integer ::IJK, nc
    nc = box%n_cell
    inv_dr = 1.0_dp/box%dr

    do j=1,nc
      do i=1,nc
        avg = 0.25*(box%cc(i+1,j,i_rho) + &
                    box%cc(i-1,j,i_rho) + &
                    box%cc(i,j+1,i_rho) + &
                    box%cc(i,j-1,i_rho))
        box%cc(i,j,i_rho) = box%cc(i,j, i_rho) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if1) - box%fc(i,j,1,if1)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if1) - box%fc(i,j,2,if1)))
       avg = 0.25*(box%cc(i+1,j,i_mom(1)) + &
                    box%cc(i-1,j,i_mom(1)) + &
                    box%cc(i,j+1,i_mom(1)) + &
                    box%cc(i,j-1,i_mom(1)))
       box%cc(i,j,i_mom(1)) = box%cc(i,j, i_mom(1)) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if2) - box%fc(i,j,1,if2)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if2) - box%fc(i,j,2,if2)))

       avg = 0.25*(box%cc(i+1,j,i_mom(2)) + &
                    box%cc(i-1,j,i_mom(2)) + &
                    box%cc(i,j+1,i_mom(2)) + &
                    box%cc(i,j-1,i_mom(2)))
       box%cc(i,j,i_mom(2)) = box%cc(i,j, i_mom(2)) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if3) - box%fc(i,j,1,if3)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if3) - box%fc(i,j,2,if3)))

       avg = 0.25*(box%cc(i+1,j,i_e) + &
                    box%cc(i-1,j,i_e) + &
                    box%cc(i,j+1,i_e) + &
                    box%cc(i,j-1,i_e))
       box%cc(i,j,i_e) = box%cc(i,j, i_e) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if4) - box%fc(i,j,1,if4)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if4) - box%fc(i,j,2,if4)))


      end do
    end do
  end subroutine updateSoln

  !======================================================
  subroutine ref_rout( box, cell_flags )
    type(box_t), intent(in) :: box
    integer, intent(out) :: cell_flags(DTIMES(box%n_cell))
    real(dp) :: diff, tol
    integer :: IJK, nc
    nc = box%n_cell
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
     end do;CLOSE_DO


  end subroutine ref_rout

  !=====================================================================
  subroutine writePrimitives( box, new_vars, n_var )
    type(box_t), intent(in):: box
    integer, intent(in) :: n_var
    real(dp) :: new_vars(DTIMES(0:box%n_cell+1),n_var)

    !XVelocity
    new_vars(DTIMES(:), 1) = box%cc(DTIMES(:), i_mom(1))/box%cc(DTIMES(:), i_rho)
    !YVelocity
    new_vars(DTIMES(:), 2) = box%cc(DTIMES(:), i_mom(2))/box%cc(DTIMES(:), i_rho)
    !Pressure
    new_vars(DTIMES(:), 3) = (Y-1.0_dp)*(box%cc(DTIMES(:), i_e) - &
                             (box%cc(DTIMES(:), i_mom(1))**2 + &
                              box%cc(DTIMES(:), i_mom(2))**2) &
                             /(2.0_dp*box%cc(DTIMES(:), i_rho)))

  end subroutine writePrimitives


end program KT_euler
