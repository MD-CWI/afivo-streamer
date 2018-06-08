#include "../afivo/src/cpp_macros_$Dd.h"
module m_init_cond_$Dd
  use m_streamer
  use m_a$D_all
  use m_geometry


  implicit none
  private

  ! Type to store initial conditions in
  type initcnd_t
     real(dp)                        :: background_density
     real(dp)                        :: stochastic_density
     integer                         :: n_cond
     real(dp), allocatable           :: seed_r0(:, :)
     real(dp), allocatable           :: seed_r1(:, :)
     real(dp), allocatable           :: seed_density(:)
     integer, allocatable            :: seed_charge_type(:)
     real(dp), allocatable           :: seed_width(:)
     character(ST_slen), allocatable :: seed_falloff(:)
     character(ST_slen)              :: die_type    
     real(dp), allocatable           :: die_center(:)
     real(dp), allocatable           :: die_axis(:)
     real(dp)                        :: geo_parameter
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected, public :: init_conds

  public :: init_cond_initialize
  public :: init_cond_set_box
  public :: init_cond_stochastic_density
  public :: define_DI
  public :: a$D_bc_dirichlet_zero_fc
  public :: a$D_loop_box_arg_DI
  public :: a$D_consistent_fluxes_DI
  public :: a$D_loop_boxes_DI
  public :: eps01_gc2
  
contains

  ! Set the initial conditions from the configuration
  subroutine init_cond_initialize(cfg, n_dim)
    type(CFG_t), intent(inout) :: cfg !< Settings
    integer, intent(in)        :: n_dim

    integer                    :: n_cond, varsize
    real(dp)                   :: dlen
    real(dp), allocatable      :: tmp_vec(:)
    type(initcnd_t)            :: ic

    call CFG_add(cfg, "background_density", 1.0e14_dp, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "stochastic_density", 0.0_dp, &
         "Stochastic background density (1/m3)")
    call CFG_add(cfg, "seed_density", [1.0e16_dp], &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed_rel_r0", [DTIMES(0.5_dp)], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", [DTIMES(0.5_dp)], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_charge_type", [0], &
         "Type of seed: neutral (0), ions (1) or electrons (-1)", .true.)
    call CFG_add(cfg, "seed_width", [0.25d-3], &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed_falloff", ['smoothstep'], &
         "Fall-off type for seed (sigmoid, gaussian, smoothstep, step, laser)", .true.)
         

    call CFG_add(cfg, "die_type", ['square'], &
         "Geometric shape of dielectric (square, elliptic, trapz)")     
    call CFG_add(cfg, "die_rel_center",[DTIMES(0.5_dp)], &
         "Relative central position of the dielectric", .true.)
    call CFG_add(cfg, "die_rel_axis",[DTIMES(0.5_dp)], &
         "Relative axis(sides or semi-axis) size of the dielectric", .true.)
    call CFG_add(cfg, "geo_parameter",[0.25_dp], &
         "Extra geometrical parameter, required for some geometries ", .true.)    




    call CFG_get_size(cfg, "seed_density", n_cond)
    ic%n_cond = n_cond

    call CFG_get_size(cfg, "seed_rel_r0", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_rel_r0 variable has incompatible size"

    call CFG_get_size(cfg, "seed_rel_r1", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_rel_r1 variable has incompatible size"

    call CFG_get_size(cfg, "seed_charge_type", varsize)
    if (varsize /= n_cond) &
         stop "seed_charge_type variable has incompatible size"

    call CFG_get_size(cfg, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_width variable has incompatible size"
         
    call CFG_get_size(cfg, "die_rel_center", varsize)
    if (varsize /= n_dim) &
         stop "die_rel_center variable has incompatible size"

    call CFG_get_size(cfg, "die_rel_axis", varsize)
    if (varsize /= n_dim) &
         stop "die_rel_axis variable has incompatible size"

    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(n_dim, n_cond))
    allocate(ic%seed_r1(n_dim, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))
    allocate(ic%die_center(n_dim))
    allocate(ic%die_axis(n_dim))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    call CFG_get(cfg, "domain_len", dlen)
    ic%seed_r0 = dlen * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = dlen * reshape(tmp_vec, [n_dim, n_cond])
   
    
    deallocate(tmp_vec)
    allocate(tmp_vec(n_dim))
    call CFG_get(cfg, "die_rel_center", tmp_vec)
    ic%die_center = dlen * tmp_vec
    call CFG_get(cfg, "die_rel_axis", tmp_vec)
    ic%die_axis = dlen * tmp_vec

    
    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "stochastic_density", ic%stochastic_density)
    call CFG_get(cfg, "seed_density", ic%seed_density)
    call CFG_get(cfg, "seed_charge_type", ic%seed_charge_type)
    call CFG_get(cfg, "seed_width", ic%seed_width)
    call CFG_get(cfg, "seed_falloff", ic%seed_falloff)
    call CFG_get(cfg, "geo_parameter", ic%geo_parameter)
    call CFG_get(cfg, "die_type", ic%die_type)
    
    init_conds = ic

  end subroutine init_cond_initialize

  !> Add a stochastic background density to the electrons and ions. Note: this
  !> routine temporarily uses variable i_rhs
  subroutine init_cond_stochastic_density(tree)
    use m_DI_prolong_$Dd
    use m_DI_restrict_$Dd
    use m_DI_ghostcell_$Dd
    type(a$D_t), intent(inout) :: tree
    integer                    :: my_lvl, lvl, i, id

    if (init_conds%stochastic_density <= 0.0_dp) return

    ! Determine at which level to create the background density. This is the
    ! highest level that is fully refined
    do my_lvl = 1, tree%highest_lvl
       if (size(tree%lvls(my_lvl)%leaves) > 0) exit
    end do

    ! Use i_rhs to store the stochastic density at this level
    call a$D_tree_clear_cc(tree, i_rhs)
    !$omp do private(id)
    do i = 1, size(tree%lvls(my_lvl)%ids)
       id = tree%lvls(my_lvl)%ids(i)
       call set_stochastic_density(tree%boxes(id))
    end do
    !$omp end do

    ! Prolong to finer levels. The coarser (hidden) levels are set at the end.
    do lvl = my_lvl, tree%highest_lvl-1
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a$D_gc_box_DI(tree%boxes, id, i_eps, i_rhs, &
               a$D_gc_interp_lim_DI, a$D_bc_neumann_zero, ST_epsilon_die)
       end do

       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a$D_prolong_linear_from_DI(tree%boxes, id, i_rhs, i_eps = i_eps, eps_DI = ST_epsilon_die, add = .true.)
       end do

    end do

    ! Finally, add the stochastic density to the electrons and ions
    do lvl = my_lvl, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a$D_box_add_cc(tree%boxes(id), i_rhs, i_electron)
          call a$D_box_add_cc(tree%boxes(id), i_rhs, i_pos_ion)
       end do
    end do

    ! Restrict and fill ghost cells
    call a$D_restrict_tree_DI(tree, i_electron, i_eps = i_eps, eps_DI = ST_epsilon_die, s_flag = .false.)
    call a$D_restrict_tree_DI(tree, i_pos_ion, i_eps = i_eps, eps_DI = ST_epsilon_die, s_flag = .false.)
    call a$D_gc_tree_DI(tree, i_eps, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero, a$D_gc_interp_lim_DI, eps_DI = ST_epsilon_die)
    call a$D_gc_tree_DI(tree, i_eps, i_pos_ion, a$D_gc_interp_lim, a$D_bc_neumann_zero, a$D_gc_interp_lim_DI, eps_DI = ST_epsilon_die)

  end subroutine init_cond_stochastic_density

  subroutine set_stochastic_density(box) 
    use omp_lib
    use m_streamer, only: ST_prng
    type(box$D_t), intent(inout) :: box
    integer                      :: proc_id, IJK
    real(dp)                     :: density, med
    
    med = a$D_harm_w(1.0_dp, ST_epsilon_die, 0.5_dp)

    proc_id = 1+omp_get_thread_num()

    do KJI_DO(1,box%n_cell)
      if(box%cc(IJK, i_eps) < ST_epsilon_die) then
        density = ST_prng%rngs(proc_id)%unif_01() * &
        init_conds%stochastic_density
        box%cc(IJK, i_rhs) = density
      else
        box%cc(IJK, i_rhs) = 0.0_dp
      end if
    end do; CLOSE_DO
  end subroutine set_stochastic_density

  !> Sets the initial condition
  subroutine init_cond_set_box(box)
    use m_geometry
    use m_units_constants

    type(box$D_t), intent(inout) :: box
    integer                     :: IJK, n, nc
    real(dp)                    :: rr($D)
    real(dp)                    :: density
    character(len=100)  :: input_file
    type(lookup_table_t)    :: FL_lkp
    type(LT_loc_t)             :: loc
    integer             :: FL_if_elec, FL_if_ion
    real(dp), allocatable   :: x_data(:), y_data(:)  

    nc = box%n_cell
    box%cc(DTIMES(:), i_electron)   = init_conds%background_density
    box%cc(DTIMES(:), i_pos_ion)    = init_conds%background_density
    box%fc(DTIMES(:), :, sigma_rhs) = 0.0_dp
    box%cc(DTIMES(:), i_phi)        = 0.0_dp ! Inital potential set to zero
    
    

    do KJI_DO(0,nc+1)
       rr   = a$D_r_cc(box, [IJK])
       do n = 1, init_conds%n_cond
         density = init_conds%seed_density(n) * &
              GM_density_line(rr, init_conds%seed_r0(:, n), &
              init_conds%seed_r1(:, n), $D, &
              init_conds%seed_width(n), &
              init_conds%seed_falloff(n))

           ! Add electrons and/or ions depending on the seed charge type
           ! (positive, negative or neutral)
           if (init_conds%seed_charge_type(n) <= 0) then
             box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + density
           end if

           if (init_conds%seed_charge_type(n) >= 0) then
             box%cc(IJK, i_pos_ion) = box%cc(IJK, i_pos_ion) + density
           end if
       end do
    end do; CLOSE_DO
    
    call define_DI(box)

  end subroutine init_cond_set_box
  
  
  subroutine define_DI(box)

    type(box$D_t), intent(inout)    :: box
    integer                         :: IJK, n, nc
    real(dp)                        :: rr($D), dr, f
    
    nc = box%n_cell
    dr = box%dr 

    do KJI_DO(0, nc+1)
       rr = a$D_r_cc(box, [IJK])
       f = vol_frac_inside(box, [IJK])
       box%cc(IJK, i_eps)        = a$D_harm_w(1.0_dp, ST_epsilon_die, f)
       if (f == 1.0_dp) then
         box%cc(IJK, i_electron)   = 0.0_dp
         box%cc(IJK, i_pos_ion)    = 0.0_dp
       end if
    end do; CLOSE_DO
 
 !   if (maxval(box%cc(DTIMES(:), i_eps)) > minval(box%cc(DTIMES(:), i_eps))) then  
      
!      do KJI_DO(0, nc+1)
!
!         if(maxval([IJK]) < nc+1 .and. minval([IJK])> 0) then
!#if $D == 2
!           if (box%cc(IJK, i_eps) == ST_epsilon_die .and. &
!              (box%cc(i-1, j, i_eps) < ST_epsilon_die .and. box%cc(i+1, j, i_eps) < ST_epsilon_die) ) then
!             box%cc(IJK, i_eps) = a$D_harm_w(box%cc(i-1, j, i_eps), box%cc(i+1, j, i_eps), 0.5_dp)
!           end if
!        
!           if (box%cc(IJK, i_eps) == ST_epsilon_die .and. &
!              (box%cc(i, j-1, i_eps) < ST_epsilon_die .and. box%cc(i, j+1, i_eps) < ST_epsilon_die) ) then
!             box%cc(IJK, i_eps) = a$D_harm_w(box%cc(i, j-1, i_eps), box%cc(i, j+1, i_eps), 0.5_dp)
!           end if
           
!#elif $D == 3

!#endif
 !        end if
  !    end do; CLOSE_DO
    
   ! end if

  end subroutine define_DI
  
  
  
  real(dp) function vol_frac_inside(box, ix)
    type(box$D_t), intent(in) :: box
    integer, intent(in)       :: ix($D)
    real(dp)                  :: dr, rr($D), r($D), in_counts
    integer                   :: IJK
    
    dr = box%dr
    rr = a$D_r_cc(box, ix)
    
    in_counts = 0.0_dp
    do KJI_DO(1, 8)
#if $D == 2
      r(1) = rr(1) + i*dr/8.0_dp - dr/2.0_dp - dr/16.0_dp
      r(2) = rr(2) + j*dr/8.0_dp - dr/2.0_dp - dr/16.0_dp
#elif $D == 3
      r(1) = rr(1) + i*dr/8.0d0 - dr/2.0d0 - dr/16.0d0
      r(2) = rr(2) + j*dr/8.0d0 - dr/2.0d0 - dr/16.0d0
      r(3) = rr(3) + k*dr/8.0d0 - dr/2.0d0 - dr/16.0d0
#endif
      if (DI_interior(init_conds%die_type, init_conds%die_center(:), &
         init_conds%die_axis(:), $D, r(:))) then
        in_counts=in_counts+1.0_dp
      end if
    end do; CLOSE_DO
    
    vol_frac_inside = in_counts / (8.0_dp)**$D
  end function vol_frac_inside
  
  function eps01_gc2(box) !> Including 2nd ghost cell, return 0 for inside, 1 for outside

    type(box$D_t), intent(in)    :: box
    integer                         :: IJK, nc
    real(dp)                        :: rr($D), dr, f 
    real(dp)                        :: eps01_gc2(DTIMES(-1:box%n_cell+2))
  
    
    nc = box%n_cell
    dr = box%dr 

    do KJI_DO(-1, nc+2)
       rr = a$D_r_cc(box, [IJK])
       f = vol_frac_inside(box, [IJK])
       if (maxval([IJK]) < nc+1 .and. minval([IJK])> 0) then
         eps01_gc2(IJK) = a$D_delta(box%cc(IJK, i_eps),ST_epsilon_die)
       else
         eps01_gc2(IJK) = a$D_delta(a$D_harm_w(1.0_dp, ST_epsilon_die, f), ST_epsilon_die)
       end if
    end do; CLOSE_DO



  end function eps01_gc2
   
  subroutine a$D_bc_dirichlet_zero_fc(box, nb, s_iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb, s_iv
    integer, intent(out)        :: bc_type

    bc_type = af_bc_dirichlet
    call a$D_set_box_gc_fc(box, nb, s_iv, gc_scalar = [DTIMES(0.0_dp)])

  end subroutine a$D_bc_dirichlet_zero_fc
  
  !> Call procedure for each box in tree (not fully inside the dielectric), with argument rarg
  subroutine a$D_loop_box_arg_DI(tree, my_procedure, rarg, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_arg)        :: my_procedure
    real(dp), intent(in)          :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             if (minval(tree%boxes(id)%cc(DTIMES(:), i_eps)) /= ST_epsilon_die) then
               call my_procedure(tree%boxes(id), rarg)
             end if
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             if (minval(tree%boxes(id)%cc(DTIMES(:), i_eps)) /= ST_epsilon_die) then
               call my_procedure(tree%boxes(id), rarg)
             end if
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_box_arg_DI
  
  subroutine a$D_consistent_fluxes_DI(tree, f_ixs)
    type(a$D_t), intent(inout)     :: tree         !< Tree to operate on
    integer, intent(in)           :: f_ixs(:)     !< Indices of the fluxes
    integer                       :: lvl, i, id, nb, nb_id

    if (.not. tree%ready) stop "Tree not ready"
    !$omp parallel private(lvl, i, id, nb, nb_id)
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          do nb = 1, a$D_num_neighbors
             nb_id = tree%boxes(id)%neighbors(nb)

             ! If the neighbor exists and has no children, set flux
             if (nb_id > af_no_box) then
                if (.not. a$D_has_children(tree%boxes(nb_id)) .and. &
                   minval(tree%boxes(nb_id)%cc(DTIMES(:), i_eps)) /= ST_epsilon_die) then
                  call flux_from_children(tree%boxes, id, nb, f_ixs)
                end if
             end if
          end do
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine a$D_consistent_fluxes_DI
  
  !> Used for flux_elec and flux_ion
  subroutine a$D_loop_boxes_DI(tree, my_procedure, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_boxes)      :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             if (minval(tree%boxes(id)%cc(DTIMES(:), i_eps)) /= ST_epsilon_die) then
               call my_procedure(tree%boxes, id)
             end if
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             if (minval(tree%boxes(id)%cc(DTIMES(:), i_eps)) /= ST_epsilon_die) then
               call my_procedure(tree%boxes, id)
             end if
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_boxes_DI
  

end module m_init_cond_$Dd