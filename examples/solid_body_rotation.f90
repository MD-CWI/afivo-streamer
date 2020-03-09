#include "../src/cpp_macros.h"

program solid_body_rotation
  use m_af_all
  implicit none
	
	integer, parameter :: ncells = 8
	integer :: iq, i_flux
	integer, parameter :: coord_type = af_xyz
	real(dp) :: l_max(NDIM), l_min(NDIM)
	integer :: grid(NDIM)
	logical :: periodicBC(NDIM)
	
	type(af_t) :: tree
	!Add the stuff for AMR later
	real(dp) :: dt, time, end_time
	integer :: time_steps, t_iter
	character(len=100) :: fname
	print *, "Running solid body rotation 2D"
	print *, "Number of threads", af_get_max_threads()
	
	!Setting up the domain and ncells
	grid(:) = 8*ncells
	l_max(:) = 1.0_dp
	l_min(:) = -1.0_dp
	periodicBC(:) = .false.
	
	
	!Initialize the variables, boundary conditions and the tree
	call af_add_cc_variable(tree, "q", ix=iq)
	call af_add_fc_variable(tree, "flux", ix=i_flux)
	call af_set_cc_methods(tree, iq, af_bc_dirichlet_zero)
	
	call af_init(tree, ncells, l_max, grid, & 
			 periodic = periodicBC, r_min = l_min, &
			 coord = coord_type)
			 
	!Setting the boundary conditions
	call af_loop_box(tree, setInitConds)
	
	!Fill the ghost cells according to the initial conditions
	call af_gc_tree(tree, [iq])
	
	!Setting the time step data
	dt = 0.005_dp
	time = 0.0
	end_time = acos(-1.0_dp)
	time_steps = ceiling(end_time/dt)
	
	!Starting the time integration
	do t_iter = 1, time_steps
		
		if (mod(t_iter, 10) == 0) then
			write(fname, "(A,I0)") "solid_body_rotation_" // DIMNAME // "_", t_iter
			call af_write_silo(tree, trim(fname), t_iter, time, dir="output")
		end if
		
		
		call af_loop_tree(tree, korenFlux)
		!print *, 'Finished Looping over tree'
		
		call af_loop_box_arg(tree, updateSoln, [dt])
		!print *, 'finished time stepping'
		
		call af_gc_tree(tree, [iq])
		!print *, 'Updating ghost cells'
		
		time = time + dt
		
		
	end do
	call af_destroy(tree)
	
	contains
	
	subroutine setInitConds( box )
		type(box_t), intent(inout) :: box
		integer :: IJK, nc
		real(dp) :: rr(NDIM)
		
		nc = box%n_cell
		do KJI_DO(0, nc+1)
			rr = af_r_cc(box, [IJK])
			box%cc(IJK, iq) = solution(rr)
		end do; CLOSE_DO
	end subroutine setInitConds
	
	function solution( rr ) result(sol)
		real(dp), intent(in) :: rr(NDIM)
		real(dp) :: sol
		real(dp) :: peak(NDIM), xLim(NDIM), yLim(NDIM)
		
		peak(1) = 0.45_dp
		peak(2) = 0.0_dp
		
		xLim(1) = 0.1_dp
		xLim(2) = 0.6_dp
		
		yLim(1) = -0.25_dp
		yLim(2) = 0.25_dp
		
		if( xLim(1) < rr(1) .and. rr(1) < xLim(2) .and. yLim(1) < rr(2) & 
		  .and. rr(2) < yLim(2) ) then
			sol = 1.0_dp
		else if (norm2(rr + peak) < 0.35_dp) then
			sol = 1.0_dp - (norm2(rr + peak)/0.35_dp)
		else
			sol = 0.0_dp
		end if	
		
	end function
	
  subroutine korenFlux( tree, id )
    use m_af_flux_schemes
    type(af_t), intent(inout)  :: tree
    integer, intent(in) :: id
    integer :: nc, IJK
    real(dp), allocatable :: cc(DTIMES(:), :)
    real(dp), allocatable :: v(DTIMES(:),:)
    real(dp) :: rr(NDIM)
		
    nc = tree%boxes(id)%n_cell
    allocate(cc(DTIMES(-1:nc+2), 1))
    allocate(v(DTIMES(1:nc+1), NDIM))
		
		!print *, "Inside Koren flux", id
    call af_gc2_box(tree, id, [iq], cc)
		!Computing the spatially varying velocities
    do KJI_DO(1,nc+1)
      rr = af_r_cc(tree%boxes(id), [IJK]) - 0.5_dp*tree%boxes(id)%dr
			v(IJK, 1) = 2.0_dp*rr(2)
      v(IJK, 2) = -2.0_dp*rr(1) 
    end do; CLOSE_DO
		
    call flux_koren_2d(cc(DTIMES(:), iq), v, nc, 2)
    tree%boxes(id)%fc(:,:,:,i_flux) = v
		
  end subroutine korenFlux
	
  subroutine updateSoln( box, dt )
    type(box_t), intent(inout) :: box
    real(dp), intent(in) :: dt(:)
    real(dp) :: inv_dr(NDIM)
    integer :: IJK, nc
		
    nc = box%n_cell
    inv_dr = 1/box%dr
		
    do j=1,nc
      do i=1,nc
        box%cc(i,j,iq) = box%cc(i,j,iq) - dt(1) * ( & 
                        inv_dr(1)* & 
                        (box%fc(i+1,j,1,iq) - box%fc(i,j,1,iq)) &
                        + inv_dr(2)* & 
                        (box%fc(i,j+1,2,iq) - box%fc(i,j,2,iq)))
      end do
    end do
		
  end subroutine updateSoln
	
end program solid_body_rotation
