#include "../src/cpp_macros.h"

program solid_body_rotation
  use m_af_all
  implicit none
	
	integer, parameter :: ncells = 8
	integer :: iq, i_flux
	integer, parameter :: coord_type = af_xyz
	real(dp) :: l_max(NDIM), l_min(NDIM), error
	integer :: grid(NDIM)
	logical :: periodicBC(NDIM)
	
	type(af_t) :: tree
	real(dp) :: dt, time, end_time
	integer :: time_steps, t_iter
	character(len=100) :: fname
	
	!AMR related stuff
	type(ref_info_t) :: refine_info
	integer :: refine_steps
	real(dp) :: dr_min(NDIM)
	
	print *, "Running solid body rotation 2D"
	print *, "Number of threads", af_get_max_threads()
	
	!Setting up the domain and ncells
	grid(:) = 4*ncells
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
			 
	!Setting the boundary conditions & fill ghost cells 
	call af_loop_box(tree, setInitConds)
	call af_gc_tree(tree, [iq])
	!One level of refinement
	call af_adjust_refinement(tree, refine_routine, refine_info, 1)
	call af_restrict_tree(tree, iq)
	call af_gc_tree(tree, [iq])
	
	!call af_write_silo(tree, "testsbr_amr", dir="output")
	!Setting the time step data
	time = 0.0
	end_time = 0.5_dp*acos(-1.0_dp)
	time_steps = ceiling(end_time/dt)
	
	!Starting the time integration
	do
	
	  dr_min = af_min_dr(tree)
	  dt = 0.5_dp/(sum(3.0_dp/dr_min) + epsilon(1.0_dp))
		
		if (mod(t_iter, 50) == 0) then
			write(fname, "(A,I0)") "amr_solid_body_rotation_" // DIMNAME // "_", t_iter
			call af_write_silo(tree, trim(fname), t_iter, time, dir="output")
		end if
		
		
		call af_loop_tree(tree, korenFlux, .true.)
		call af_consistent_fluxes(tree, [iq])
		!print *, 'Finished Looping over tree'
		
		call af_loop_box_arg(tree, updateSoln, [dt])
		!print *, 'finished time stepping'
		call af_restrict_tree(tree, iq)
		
		call af_gc_tree(tree, [iq])
		!print *, 'Updating ghost cells'
		
		call af_adjust_refinement(tree, refine_routine, refine_info, 1)
		
		time = time + dt
		t_iter = t_iter+1
		
		!Checking for divergence
		call af_tree_maxabs_cc(tree, iq, error)
		if (error > 5.0_dp .or. dt < 1e-6_dp) then
		  print *, "Solution diverging!"
		  exit
		 end if
		
		if (time > end_time) exit
		
	end do
	call af_destroy(tree)
	!==============================================================================================================================
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
		
	end function solution
	
	!=============================================================================================
	subroutine refine_routine( box, cell_flags )
	  type(box_t), intent(in) :: box
	  integer, intent(out) :: cell_flags(DTIMES(box%n_cell))
	  real(dp) :: diff
	  integer :: IJK, nc
	  
	  nc = box%n_cell
	  do KJI_DO(1,nc)
      diff = box%dr(1)**2*abs(box%cc(i+1,j,iq)+box%cc(i-1,j,iq)-2*box%cc(i,j,iq)) + &
              box%dr(2)**2*abs(box%cc(i,j+1,iq)+box%cc(i,j-1,iq)-2*box%cc(i,j,iq))
      if (diff > 1.0e-5_dp .and. box%lvl .le. 3) then
        cell_flags(IJK) = af_do_ref
      else if (diff < 0.1_dp*1.0e-5_dp) then
        cell_flags(IJK) = af_rm_ref
      else
        cell_flags(IJK) = af_keep_ref
      end if	    
	  end do; CLOSE_DO
	end subroutine refine_routine
	!=============================================================================================
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
		
    call flux_koren_2d(cc(DTIMES(:), 1), v, nc, 2)
    tree%boxes(id)%fc(:,:,:,i_flux) = v
		
  end subroutine korenFlux
	!=============================================================================================	
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
