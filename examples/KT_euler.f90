#include "../src/cpp_macros.h"

program KT_euler
  use m_af_all
  implicit none
  
  integer, parameter :: ncells = 8
  integer :: ic1, ic2, ic3, ic4
  integer :: if1, if2, if3, if4
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
  
  
  
  grid(:) = 25*ncells
  l_max(:) = 1.0_dp
  l_min(:) = 0.0_dp
  periodicBC(:) = .false.
  
  call af_add_cc_variable(tree, "rho", ix=ic1)
  call af_add_cc_variable(tree, "rhoU", ix=ic2)
  call af_add_cc_variable(tree, "rhoV", ix=ic3)
  call af_add_cc_variable(tree, "E", ix=ic4)
  
  call af_add_fc_variable(tree, "flux1", ix=if1)
  call af_add_fc_variable(tree, "flux2", ix=if2)
  call af_add_fc_variable(tree, "flux3", ix=if3)
  call af_add_fc_variable(tree, "flux4", ix=if4)
  
  call af_set_cc_methods(tree, ic1, af_bc_neumann_zero)
  call af_set_cc_methods(tree, ic2, af_bc_neumann_zero)
  call af_set_cc_methods(tree, ic3, af_bc_neumann_zero)
  call af_set_cc_methods(tree, ic4, af_bc_neumann_zero)
  
  call af_init(tree, ncells, l_max, grid, &
               periodic = periodicBC, r_min= l_min, &
               coord = coord_type)

  
  !Init mesh refinement
!  do 
!      refine_steps = refine_steps + 1
!      !Settng init conds for each refinement is needed as we use that data as
!      !refinement criterion
!      call af_loop_box(tree, setInitConds)
!      call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
!      call af_adjust_refinement(tree, ref_rout, refine_info, 1)
!  if (refine_info%n_add == 0) exit
!  end do
!  call af_restrict_tree(tree, ic1)
!  call af_restrict_tree(tree, ic2)
!  call af_restrict_tree(tree, ic3)
!  call af_restrict_tree(tree, ic4)
!  call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
  
               
  call af_loop_box(tree, setInitConds)
  
  call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
  
  !call af_write_silo(tree, 'eulerInit', dir='output')
  
  !Setting the timestep data
  time = 0.0_dp
  end_time = 0.2_dp
  t_iter = 0
  do 
   !Updating the primitive vars
   !call af_tree_copy_cc(tree, ic2, ip2)
   !call af_tree_copy_cc(tree, ic3, ip3)
   !call af_tree_copy_cc(tree, ic4, ip4)
   !call af_loop_box(tree, updatePrimitives)
   dr_min = af_min_dr(tree)
   dt = end_time/(sum(wSp/dr_min) + epsilon(1.0_dp))
    if (mod(t_iter, 10) == 0) then
      write(fname, "(A,I0)") "KT_euler_" // DIMNAME // "_", t_iter
      call af_write_silo(tree, trim(fname), t_iter, time, dir="output", &
           add_vars = writePrimitives, add_names=["xVel","yVel","pres"])
    end if
  
   call af_loop_tree(tree, fluxComputation)
   call af_consistent_fluxes(tree, [ic1,ic2, ic3, ic4])
   call af_loop_box_arg(tree, updateSoln, [dt])
   call af_restrict_tree(tree, ic1)
   call af_restrict_tree(tree, ic2)
   call af_restrict_tree(tree, ic3)
   call af_restrict_tree(tree, ic4)
   call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
   !call af_adjust_refinement(tree, ref_rout, refine_info, 1)
   
   !call af_loop_box_arg(tree, updateSoln, [dt])
   !call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
    
   t_iter = t_iter + 1
   time = time + dt
   
   
   call af_tree_maxabs_cc(tree, ic1, error)
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
  !======================================================
!  subroutine updatePrimitives( box )
!    type(box_t), intent(inout) :: box
!    real(dp) :: rho, rhou, rhov, E
!    integer :: nc, IJK
!    nc = box%n_cell
!    do KJI_DO(1,nc)
!      rho = box%cc(IJK, ic1)
!      rhou = box%cc(IJK, ic2)
!      rhov = box%cc(IJK, ic3)
!      E = box%cc(IJK, ic4)
!      box%cc(IJK, ip2) = rhou/rho
!      box%cc(IJK, ip3) = rhov/rho
!      box%cc(IJK, ip4) = (Y-1.0_dp)*(E - & 
!                                    (rhou**2 + rhov**2)/(2.0_dp*rho))  
!    end do; CLOSE_DO

!    !box%cc(DTIMES(:),ip2) = box%cc(DTIMES(:),ip2)/box%cc(DTIMES(:),ic1)
!    !box%cc(DTIMES(:),ip2) = box%cc(DTIMES(:),ip3)/box%cc(DTIMES(:),ic1)
!     
!    
!    
!    
!  end subroutine updatePrimitives
  !======================================================
!  subroutine updatePrimitives( tree, id )
!    type(af_t), intent(inout) :: tree
!    integer, intent(in) :: id
!    real(dp) :: rho, rhou, rhov, E
!    integer :: nc, IJK
!    
!    nc = tree%boxes(id)%n_cell
!    do KJI_DO(1,nc)
!      rho = tree%boxes(id)%cc(IJK, ic1)
!      rhou = tree%boxes(id)%cc(IJK, ic2)
!      rhov = tree%boxes(id)%cc(IJK, ic3)
!      E = tree%boxes(id)%cc(IJK, ic4)
!      tree%boxes(id)%cc(IJK, ip2) = rhou/rho
!      tree%boxes(id)%cc(IJK, ip3) = rhov/rho
!      tree%boxes(id)%cc(IJK, ip4) = (Y-1.0_dp)*(E - & 
!                                    (rhou**2 + rhov**2)/(2.0_dp*rho))  
!    end do; CLOSE_DO
!    
!    
!    
!  end subroutine updatePrimitives
  !=====================================================================
  
  subroutine setInitConds( box )
    type(box_t), intent(inout) :: box
    integer :: IJK, nc
    real(dp) :: rr(NDIM)
    !real(dp) :: p(4), rho(4), u(4), v(4)
    real(dp) :: conservatives(4)
    
    
    nc = box%n_cell
    do KJI_DO(0, nc+1)
      rr = af_r_cc(box, [IJK])
      if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
        conservatives = convert_to_conservatives(rho(1), u(1), v(1), p(1))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
!        box%cc(IJK, ip2) = u(1)
!        box%cc(IJK, ip3) = v(1)
!        box%cc(IJK, ip4) = p(1)
      elseif (rr(1) .le. 0.5_dp .and. rr(2) .ge. 0.5_dp) then
        conservatives = convert_to_conservatives(rho(2), u(2), v(2), p(2))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
!        box%cc(IJK, ip2) = u(2)
!        box%cc(IJK, ip3) = v(2)
!        box%cc(IJK, ip4) = p(2)
      elseif (rr(1) .le. 0.5_dp .and. rr(2) .le. 0.5_dp) then
        conservatives = convert_to_conservatives(rho(3), u(3), v(3), p(3))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
!        box%cc(IJK, ip2) = u(3)
!        box%cc(IJK, ip3) = v(3)
!        box%cc(IJK, ip4) = p(3)
      else
        conservatives = convert_to_conservatives(rho(4), u(4), v(4), p(4))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
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
    integer, intent(in) :: id
    integer :: nc
    real(dp) :: dr(NDIM)
    real(dp), allocatable :: cc(DTIMES(:), :)
    real(dp), allocatable :: fc1(DTIMES(:), :), &
                             fc2(DTIMES(:), :), &
                             fc3(DTIMES(:), :), &
                             fc4(DTIMES(:), :)
    real(dp), allocatable :: locWSp(DTIMES(:), :)
    
    nc = tree%boxes(id)%n_cell
    dr = tree%boxes(id)%dr !Not used
    allocate(cc(DTIMES(-1:nc+2), 4))
    allocate(fc1(DTIMES(1:nc+1), 2*NDIM))
    allocate(fc2(DTIMES(1:nc+1), 2*NDIM))
    allocate(fc3(DTIMES(1:nc+1), 2*NDIM))
    allocate(fc4(DTIMES(1:nc+1), 2*NDIM))
    
    allocate(locWSp(DTIMES(1:nc+1), NDIM))
    !cc = 0
    call af_gc2_box(tree, id, [ic1, ic2, ic3, ic4], cc)
    fc1 = 1.0_dp
    fc2 = 1.0_dp
    fc3 = 1.0_dp
    fc4 = 1.0_dp
    
    !Computing the left, right, bottom and top cell face values
    call flux_kurganovTadmor_2d(cc(DTIMES(:), ic1), fc1, nc, 2)
    call flux_kurganovTadmor_2d(cc(DTIMES(:), ic2), fc2, nc, 2)
    call flux_kurganovTadmor_2d(cc(DTIMES(:), ic3), fc3, nc, 2)
    call flux_kurganovTadmor_2d(cc(DTIMES(:), ic4), fc4, nc, 2)
    
    !print *, fc1(1,:,1)
    !print *, "new box"
    !Computing the wave speed at the cell faces
    
    call waveSpeed_2d(fc1, fc2, fc3, fc4, locWSp, nc)
    !locWSp = 10.0_dp
    !if (id == 1) then
    !print *, 'ID:', id, "Local wave speed", locWSp
    !if (maxval(locWsp) == INFINITY) then
    !print *, "Maxvalue of speed", maxval(locWSp)
    !end if 
    !end if    
    
    
    tree%boxes(id)%fc(:,:,1,if1) = avgFluxLR(fc1,fc2,fc3,fc4,11, nc) - & 
                                   0.5_dp*locWSp(DTIMES(:),1)* &
                                   (fc1(DTIMES(:), 2) - fc1(DTIMES(:), 1))
                                   
    tree%boxes(id)%fc(:,:,2,if1) = avgFluxTB(fc1,fc2,fc3,fc4,12, nc) - & 
                                   0.5_dp*locWSp(DTIMES(:),2)* &
                                   (fc1(DTIMES(:), 4) - fc1(DTIMES(:), 3))
    
    
    tree%boxes(id)%fc(:,:,1,if2) = avgFluxLR(fc1,fc2,fc3,fc4,21, nc) - &
                                   0.5_dp*locWSp(DTIMES(:),1)* & 
                                   (fc2(DTIMES(:), 2) - fc2(DTIMES(:), 1))
    tree%boxes(id)%fc(:,:,2,if2) = avgFluxTB(fc1,fc2,fc3,fc4,22, nc) - &
                                   0.5_dp*locWSp(DTIMES(:),2)* & 
                                   (fc2(DTIMES(:), 4) - fc2(DTIMES(:), 3))
    
    
    tree%boxes(id)%fc(:,:,1,if3) = avgFluxLR(fc1,fc2,fc3,fc4,31, nc) - &
                                   0.5_dp*locWSp(DTIMES(:),1)* & 
                                   (fc3(DTIMES(:), 2) - fc3(DTIMES(:), 1))
    tree%boxes(id)%fc(:,:,2,if3) = avgFluxTB(fc1,fc2,fc3,fc4,32, nc) - &
                                   0.5_dp*locWSp(DTIMES(:),2)* & 
                                   (fc3(DTIMES(:), 4) - fc3(DTIMES(:), 3))
                                    
                                    
    tree%boxes(id)%fc(:,:,1,if4) = avgFluxLR(fc1,fc2,fc3,fc4,41, nc) - &
                                   0.5_dp*locWSp(DTIMES(:),1)* & 
                                   (fc4(DTIMES(:), 2) - fc4(DTIMES(:), 1))
                                     
    tree%boxes(id)%fc(:,:,2,if4) = avgFluxTB(fc1,fc2,fc3,fc4,42, nc) - &
                                   0.5_dp*locWSp(DTIMES(:),2)* & 
                                   (fc4(DTIMES(:), 4) - fc4(DTIMES(:), 3))

    
  end subroutine fluxComputation
  !=====================================================================
  function avgFluxLR(fc1, fc2, fc3, fc4, eqno, nc) result(avgFlux)
    integer, intent(in) :: nc
    real(dp), intent(in), dimension(DTIMES(1:nc+1), 2*NDIM) :: fc1, fc2, fc3, fc4
    integer, intent(in) :: eqno
    real(dp) :: avgFlux(DTIMES(1:nc+1))
    
    avgFlux = 0.5_dp*( eulerFlux(fc1(DTIMES(:),1), fc2(DTIMES(:),1), & 
                                 fc3(DTIMES(:),1), fc4(DTIMES(:),1), eqno) + &
                       eulerFlux(fc1(DTIMES(:),2), fc2(DTIMES(:),2), & 
                                 fc3(DTIMES(:),2), fc4(DTIMES(:),2), eqno))
  end function avgFluxLR
  
  !=====================================================================
  function avgFluxTB(fc1, fc2, fc3, fc4, eqno, nc) result(avgFlux)
    integer, intent(in) :: nc
    real(dp), intent(in), dimension(DTIMES(1:nc+1), 2*NDIM) :: fc1, fc2, fc3, fc4
    integer, intent(in) :: eqno
    real(dp) :: avgFlux(DTIMES(1:nc+1))
    
    avgFlux = 0.5_dp*( eulerFlux(fc1(DTIMES(:),3), fc2(DTIMES(:),3), & 
                                 fc3(DTIMES(:),3), fc4(DTIMES(:),3), eqno) + &
                       eulerFlux(fc1(DTIMES(:),4), fc2(DTIMES(:),4), & 
                                 fc3(DTIMES(:),4), fc4(DTIMES(:),4), eqno))
                                 
  end function avgFluxTB
  !=====================================================================
  elemental function eulerFlux(c1, c2, c3, c4, eqno) result(flux)
    real(dp), intent(in) :: c1, c2, c3, c4
    integer, intent(in) :: eqno
    real(dp) :: flux
    
    select case( eqno )
      case( 11 )
        flux = c2
      case( 12 )
        flux = c3
      case( 21 )
        flux = (c2**2/c1) + (Y-1.0_dp)*(c4 - (c2**2 + c3**2)/(2.0_dp*c1))
      case( 22 )
        flux = (c2*c3)/c1
      case( 31 )
        flux = (c2*c3)/c1
      case( 32 )
        flux = (c3**2/c1) + (Y-1.0_dp)*(c4 - (c2**2 + c3**2)/(2.0_dp*c1))
      case( 41 )
        flux = (c2/c1)*(c4 + (Y-1.0_dp)*(c4 - (c2**2 + c3**2)/(2.0_dp*c1)))
      case( 42 )
        flux = (c3/c1)*(c4 + (Y-1.0_dp)*(c4 - (c2**2 + c3**2)/(2.0_dp*c1)))
    end select
  
  end function eulerFlux
  
  !=====================================================================
  subroutine waveSpeed_2d( fc1, fc2, fc3, fc4, a, nc )
    integer, intent(in) :: nc
    real(dp), intent(in) :: fc1(DTIMES(1:nc+1), 4), &
                            fc2(DTIMES(1:nc+1), 4), & 
                            fc3(DTIMES(1:nc+1), 4), & 
                            fc4(DTIMES(1:nc+1), 4) 
    real(dp), intent(inout) :: a(DTIMES(1:nc+1), 2)
    real(dp) :: u_int(1:nc+1), u_out(1:nc+1), v_int(1:nc+1), v_out(nc+1)
    real(dp) :: p_int(1:nc+1), p_out(1:nc+1), c_int(1:nc+1), c_out(1:nc+1)
    integer :: n
   
    do n=1,nc+1
      !x-dir
      !print *, product(fc1(:,n,1))
      u_int(:) = fc2(:,n,1)/fc1(:,n,1)
      u_out(:) = fc2(:,n,2)/fc1(:,n,2)
      p_int(:) = (Y-1.0_dp)*(fc4(:,n,1) - (fc2(:,n,1)**2 + fc3(:,n,1)**2)/ & 
                                          (2.0_dp*fc1(:,n,1)))
      p_out(:) = (Y-1.0_dp)*(fc4(:,n,2) - (fc2(:,n,2)**2 + fc3(:,n,2)**2)/ & 
                                          (2.0_dp*fc1(:,n,2)))
      c_int(:) = sqrt((Y * p_int)/fc1(:,n,1))
      c_out(:) = sqrt((Y * p_out)/fc1(:,n,2))
      a(:,n, 1) = max( & 
                  max(abs(u_int(:) + c_int(:)), abs(u_int(:) - c_int(:)), abs(c_int(:))), &
                  max(abs(u_out(:) + c_out(:)), abs(u_out(:) - c_out(:)), abs(c_out(:))))
                  
                  
      v_int(:) = fc2(n,:,1)/fc1(n,:,1)
      v_out(:) = fc2(n,:,2)/fc1(n,:,2)
      p_int(:) = (Y-1.0_dp)*(fc4(n,:,1) - (fc2(n,:,1)**2 + fc3(n,:,1)**2)/ & 
                                          (2.0_dp*fc1(:,n,1)))
      p_out(:) = (Y-1.0_dp)*(fc4(n,:,2) - (fc2(n,:,2)**2 + fc3(n,:,2)**2)/ & 
                                          (2.0_dp*fc1(n,:,2)))
      c_int(:) = sqrt((Y - p_int)/fc1(n,:,1))
      c_out(:) = sqrt((Y - p_out)/fc1(n,:,2))
      a(:,n, 2) = max( & 
                  max(abs(v_int(:) + c_int(:)), abs(v_int(:) - c_int(:)), abs(c_int(:))), &
                  max(abs(v_out(:) + c_out(:)), abs(v_out(:) - c_out(:)), abs(c_out(:))))
      
    end do
    
  
  end subroutine waveSpeed_2d
  
  !=====================================================================
  
  subroutine updateSoln( box, dt )
    type(box_t), intent(inout) :: box
    real(dp), intent(in) :: dt(:)
    real(dp) :: inv_dr(NDIM), avg
    integer ::IJK, nc
    nc = box%n_cell
    inv_dr = 1.0_dp/box%dr
  
    do j=1,nc
      do i=1,nc
        avg = 0.25*(box%cc(i+1,j,ic1) + &
                    box%cc(i-1,j,ic1) + &
                    box%cc(i,j+1,ic1) + &
                    box%cc(i,j-1,ic1))
        box%cc(i,j,ic1) = box%cc(i,j, ic1) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if1) - box%fc(i,j,1,if1)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if1) - box%fc(i,j,2,if1)))
       avg = 0.25*(box%cc(i+1,j,ic2) + &
                    box%cc(i-1,j,ic2) + &
                    box%cc(i,j+1,ic2) + &
                    box%cc(i,j-1,ic2))
       box%cc(i,j,ic2) = box%cc(i,j, ic2) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if2) - box%fc(i,j,1,if2)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if2) - box%fc(i,j,2,if2)))
      
       avg = 0.25*(box%cc(i+1,j,ic3) + &
                    box%cc(i-1,j,ic3) + &
                    box%cc(i,j+1,ic3) + &
                    box%cc(i,j-1,ic3))
       box%cc(i,j,ic3) = box%cc(i,j, ic3) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if3) - box%fc(i,j,1,if3)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if3) - box%fc(i,j,2,if3)))
       
       avg = 0.25*(box%cc(i+1,j,ic4) + &
                    box%cc(i-1,j,ic4) + &
                    box%cc(i,j+1,ic4) + &
                    box%cc(i,j-1,ic4))
       box%cc(i,j,ic4) = box%cc(i,j, ic4) - dt(1) *( &
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
      diff =  box%dr(1)**2*abs(box%cc(i+1,j,ic1)+box%cc(i-1,j,ic1) & 
                               -2*box%cc(i,j,ic1)) + &
                         box%dr(2)**2*abs(box%cc(i,j+1,ic1)+box%cc(i,j-1,ic1) &
                               -2*box%cc(i,j,ic1))
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
    new_vars(DTIMES(:), 1) = box%cc(DTIMES(:), ic2)/box%cc(DTIMES(:), ic1)
    !YVelocity
    new_vars(DTIMES(:), 2) = box%cc(DTIMES(:), ic3)/box%cc(DTIMES(:), ic1)
    !Pressure
    new_vars(DTIMES(:), 3) = (Y-1.0_dp)*(box%cc(DTIMES(:), ic4) - &
                             (box%cc(DTIMES(:), ic2)**2 + & 
                              box%cc(DTIMES(:), ic3)**2) & 
                             /(2.0_dp*box%cc(DTIMES(:), ic1)))
    
  end subroutine writePrimitives
  

end program KT_euler
