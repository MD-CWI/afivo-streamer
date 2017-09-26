module m_init_cond_2d
  use m_streamer
  use m_a2_all
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
     integer                         :: n_die
     real(dp), allocatable           :: die_density(:)
     real(dp), allocatable           :: die_center(:,:)
     real(dp), allocatable           :: die_axis(:,:)
     character(ST_slen), allocatable :: die_type(:)
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected, public :: init_conds

  public :: init_cond_initialize
  public :: init_cond_set_box
  public :: init_cond_stochastic_density
  public :: define_DI
  public :: a2_bc_dirichlet_mirror
  public :: a2_bc_dirichlet_zero_fc
  public :: a2_bc2_dirichlet_mirror
  
contains

  ! Set the initial conditions from the configuration
  subroutine init_cond_initialize(cfg, n_dim)
    type(CFG_t), intent(inout) :: cfg !< Settings
    integer, intent(in)        :: n_dim

    integer                    :: n_cond, n_die, varsize
    real(dp)                   :: dlen
    real(dp), allocatable      :: tmp_vec(:)
    type(initcnd_t)            :: ic

    call CFG_add(cfg, "background_density", 1.0e14_dp, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "stochastic_density", 0.0_dp, &
         "Stochastic background density (1/m3)")
    call CFG_add(cfg, "seed_density", [1.0e16_dp], &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed_rel_r0", [0.5_dp, 0.5_dp], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", [0.5_dp, 0.5_dp], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_charge_type", [0], &
         "Type of seed: neutral (0), ions (1) or electrons (-1)", .true.)
    call CFG_add(cfg, "seed_width", [0.25d-3], &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed_falloff", ['smoothstep'], &
         "Fall-off type for seed (sigmoid, gaussian, smoothstep, step, laser)", .true.)
         
    call CFG_add(cfg, "die_rel_center",[0.5_dp, 0.5_dp], &
         "Relative central position of the dielectric", .true.)
    call CFG_add(cfg, "die_rel_axis",[0.5_dp, 0.5_dp], &
         "Relative axis(sides or semi-axis) size of the dielectric", .true.)
    call CFG_add(cfg, "die_type", ['elliptic'], &
         "Geometric shape of dielectric (square, elliptic)", .true.)
    call CFG_add(cfg, "die_density", [0.0e1_dp], &
         "Density of the dielectric(1/m3)", .true.)

    call CFG_get_size(cfg, "seed_density", n_cond)
    ic%n_cond = n_cond
    call CFG_get_size(cfg, "die_density", n_die)
    ic%n_die = n_die

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
    if (varsize /= n_dim * n_die) &
         stop "die_rel_center variable has incompatible size"

    call CFG_get_size(cfg, "die_rel_axis", varsize)
    if (varsize /= n_dim * n_die) &
         stop "die_rel_axis variable has incompatible size"

    allocate(ic%seed_density(n_cond))
    allocate(ic%die_density(n_die))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(n_dim, n_cond))
    allocate(ic%seed_r1(n_dim, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))
    allocate(ic%die_type(n_die))
    allocate(ic%die_center(n_dim, n_die))
    allocate(ic%die_axis(n_dim, n_die))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    call CFG_get(cfg, "domain_len", dlen)
    ic%seed_r0 = dlen * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = dlen * reshape(tmp_vec, [n_dim, n_cond])
   
    deallocate(tmp_vec)
    allocate(tmp_vec(n_dim * n_die))
    call CFG_get(cfg, "die_rel_center", tmp_vec)
    ic%die_center = dlen * reshape(tmp_vec, [n_dim, n_die])
    call CFG_get(cfg, "die_rel_axis", tmp_vec)
    ic%die_axis = dlen * reshape(tmp_vec, [n_dim, n_die])

    
    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "stochastic_density", ic%stochastic_density)
    call CFG_get(cfg, "seed_density", ic%seed_density)
    call CFG_get(cfg, "die_density", ic%die_density)
    call CFG_get(cfg, "seed_charge_type", ic%seed_charge_type)
    call CFG_get(cfg, "seed_width", ic%seed_width)
    call CFG_get(cfg, "seed_falloff", ic%seed_falloff)
    call CFG_get(cfg, "die_type", ic%die_type)

    init_conds = ic

  end subroutine init_cond_initialize

  !> Add a stochastic background density to the electrons and ions. Note: this
  !> routine temporarily uses variable i_rhs
  subroutine init_cond_stochastic_density(tree)
    use m_a2_ghostcell, only: a2_bc_neumann_zero
    use m_a2_prolong
    type(a2_t), intent(inout) :: tree
    integer                    :: my_lvl, lvl, i, id

    if (init_conds%stochastic_density <= 0.0_dp) return

    ! Determine at which level to create the background density. This is the
    ! highest level that is fully refined
    do my_lvl = 1, tree%highest_lvl
       if (size(tree%lvls(my_lvl)%leaves) > 0) exit
    end do

    ! Use i_rhs to store the stochastic density at this level
    call a2_tree_clear_cc(tree, i_rhs)
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
          call a2_gc_box(tree%boxes, id, i_eps, i_rhs, &
               a2_gc_interp, a2_bc_neumann_zero)
       end do

       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_prolong_linear_from(tree%boxes, id, i_eps, i_rhs, add=.true.)
       end do

    end do

    ! Finally, add the stochastic density to the electrons and ions
    do lvl = my_lvl, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a2_box_add_cc(tree%boxes(id), i_rhs, i_electron)
          call a2_box_add_cc(tree%boxes(id), i_rhs, i_pos_ion)
       end do
    end do

    ! Restrict and fill ghost cells
    call a2_restrict_tree(tree, i_electron)
    call a2_restrict_tree(tree, i_pos_ion)
    call a2_gc_tree(tree, i_eps, i_electron, a2_gc_interp_lim, a2_bc_neumann_zero)
    call a2_gc_tree(tree, i_eps, i_pos_ion, a2_gc_interp_lim, a2_bc_neumann_zero)

  end subroutine init_cond_stochastic_density

  subroutine set_stochastic_density(box) 
    use omp_lib
    use m_streamer, only: ST_prng
    type(box2_t), intent(inout) :: box
    integer                      :: proc_id, i, j
    real(dp)                     :: density, med
    
    med = a2_harm_w(1.0_dp, ST_epsilon_die, 0.5_dp)

    proc_id = 1+omp_get_thread_num()

    do j = 1, box%n_cell; do i = 1, box%n_cell
      if(box%cc(i, j, i_eps) <= med) then
        density = ST_prng%rngs(proc_id)%unif_01() * &
        init_conds%stochastic_density
        box%cc(i, j, i_rhs) = density
      else
        box%cc(i, j, i_rhs) = 0.0_dp
      end if
    end do; end do
  end subroutine set_stochastic_density

  !> Sets the initial condition
  subroutine init_cond_set_box(box)
    use m_geometry
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, n, nc
    real(dp)                    :: rr(2)
    real(dp)                    :: density

    nc = box%n_cell
    box%cc(:, :, i_electron)   = init_conds%background_density
    box%cc(:, :, i_pos_ion)    = init_conds%background_density
    box%fc(:, :, :, sigma_rhs) = 0.0_dp
    box%cc(:, :, i_phi)        = 0.0_dp ! Inital potential set to zero
    box%cc(:, :, i_eps)        = 1.0_dp ! Base permitivity

    do j = 0, nc+1; do i = 0, nc+1
       rr   = a2_r_cc(box, [i, j])

       do n = 1, init_conds%n_cond
          density = init_conds%seed_density(n) * &
               GM_density_line(rr, init_conds%seed_r0(:, n), &
               init_conds%seed_r1(:, n), 2, &
               init_conds%seed_width(n), &
               init_conds%seed_falloff(n))

          ! Add electrons and/or ions depending on the seed charge type
          ! (positive, negative or neutral)
          if (init_conds%seed_charge_type(n) <= 0) then
             box%cc(i, j, i_electron) = box%cc(i, j, i_electron) + density
          end if

          if (init_conds%seed_charge_type(n) >= 0) then
             box%cc(i, j, i_pos_ion) = box%cc(i, j, i_pos_ion) + density
          end if
       end do
    end do; end do
    
    call define_DI(box)

  end subroutine init_cond_set_box
  
  
  subroutine define_DI(box)

    type(box2_t), intent(inout)    :: box
    integer                         :: i, j, n, nc
    real(dp)                        :: rr(2), dr, f, med
    
    nc = box%n_cell
    dr = box%dr 
    box%cc(:, :, i_eps)      = 1.0_dp ! Base permitivity
    med = a2_harm_w(1.0_dp, ST_epsilon_die, 0.5_dp)

    do j = 0,  nc+1; do i = 0,  nc+1
       rr = a2_r_cc(box, [i, j])
       f = vol_frac_inside(box, [i, j])
       do n = 1, init_conds%n_die
         if(DI_interior(init_conds%die_type(n), init_conds%die_center(:,n), &
            init_conds%die_axis(:,n), 2,rr(:))) then
               box%cc(i, j, i_eps)        = a2_harm_w(1.0_dp, ST_epsilon_die, f)
               box%cc(i, j, i_electron)   = init_conds%die_density(n)
               box%cc(i, j, i_pos_ion)    = init_conds%die_density(n)
         end if
       end do
    end do; end do
 
    if (maxval(box%cc(:, :, i_eps)) > minval(box%cc(:, :, i_eps))) then  
      
      do j = 0,  nc+1; do i = 0,  nc+1

         if(maxval([i, j]) < nc+1 .and. minval([i, j])> 0) then

           if(box%cc(i, j, i_eps) > med .and. (box%cc(i-1, j, i_eps) <= med .and. box%cc(i+1, j, i_eps) <= med) ) then
             box%cc(i, j, i_eps) = a2_harm_w(box%cc(i-1, j, i_eps), box%cc(i+1, j, i_eps), 0.5_dp)
           end if
        
           if(box%cc(i, j, i_eps) > med .and. (box%cc(i, j-1, i_eps) <= med .and. box%cc(i, j+1, i_eps) <= med) ) then
             box%cc(i, j, i_eps) = a2_harm_w(box%cc(i, j-1, i_eps), box%cc(i, j+1, i_eps), 0.5_dp)
           end if
           

         end if
      end do; end do
    
    end if

  end subroutine define_DI
  
  
  
  real(dp) function vol_frac_inside(box, ix)
    type(box2_t), intent(in) :: box
    integer, intent(in)       :: ix(2)
    real(dp)                  :: dr, rr(2), r(2), in_counts
    integer                   :: i, j, n
    
    dr = box%dr
    rr = a2_r_cc(box, ix)
    
    in_counts = 0.0_dp
    do j = 1,  8; do i = 1,  8

      r(1) = rr(1) + i*dr/8.0_dp - dr/2.0_dp - dr/16.0_dp
      r(2) = rr(2) + j*dr/8.0_dp - dr/2.0_dp - dr/16.0_dp

        do n = 1, init_conds%n_die
          if (DI_interior(init_conds%die_type(n), init_conds%die_center(:,n), &
            init_conds%die_axis(:,n), 2, r(:))) then
              in_counts=in_counts+1.0_dp
          end if
        end do
    end do; end do
    
    vol_frac_inside = in_counts / (8.0_dp)**2
  end function vol_frac_inside
  
  
  
  

  subroutine a2_bc_dirichlet_mirror(box, nb, iv, bc_type, i_eps)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb, iv
    integer, intent(out)        :: bc_type
    integer, intent(in), optional   :: i_eps

    bc_type = af_bc_dirichlet
    call a2_set_box_gc(box, nb, iv, init_conds%background_density)

  end subroutine a2_bc_dirichlet_mirror
  
  subroutine a2_bc2_dirichlet_mirror(boxes, id, nb, i_eps, iv, gc_side, nc, med)
    type(box2_t), intent(inout)   :: boxes(:)
    integer, intent(in)            :: id, nb, iv, nc
    integer, intent(in), optional  :: i_eps
    real(dp), intent(in), optional :: med

    real(dp), intent(out)          :: gc_side(nc)

    select case (nb)

    case (a2_neighb_lowx)
       gc_side = -boxes(id)%cc(2, 1:nc, iv) + 2.0_dp * init_conds%background_density
    case (a2_neighb_highx)
       gc_side = -boxes(id)%cc(nc-1, 1:nc, iv) + 2.0_dp * init_conds%background_density
    case (a2_neighb_lowy)
       gc_side = -boxes(id)%cc(1:nc, 2, iv) + 2.0_dp * init_conds%background_density
    case (a2_neighb_highy)
       gc_side = -boxes(id)%cc(1:nc, nc-1, iv) + 2.0_dp * init_conds%background_density

    end select

  end subroutine a2_bc2_dirichlet_mirror
  
  subroutine a2_bc_dirichlet_zero_fc(box, nb, s_iv, bc_type, i_eps)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb, s_iv
    integer, intent(out)        :: bc_type
    integer, intent(in), optional   :: i_eps

    bc_type = af_bc_dirichlet
    call a2_set_box_gc_fc(box, nb, s_iv, gc_scalar = [0.0_dp, 0.0_dp])

  end subroutine a2_bc_dirichlet_zero_fc
  
  

end module m_init_cond_2d

