#include "../afivo/src/cpp_macros_$Dd.h"
module m_init_cond_$Dd
  use m_streamer
  use m_a$D_all

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
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected, public :: init_conds

  public :: init_cond_initialize
  public :: init_cond_set_box
  public :: init_cond_stochastic_density

contains

  ! Set the initial conditions from the configuration
  subroutine init_cond_initialize(cfg, n_dim)
    type(CFG_t), intent(inout) :: cfg !< Settings
    integer, intent(in)        :: n_dim

    integer                    :: n_cond, varsize
    real(dp)                   :: dlen
    real(dp), allocatable      :: tmp_vec(:)
    type(initcnd_t)            :: ic

    call CFG_add(cfg, "background_density", 0.0_dp, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "stochastic_density", 0.0_dp, &
         "Stochastic background density (1/m3)")
    call CFG_add(cfg, "seed_density", [1.0e15_dp], &
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
         "Fall-off type for seed (sigmoid, gaussian, smoothstep, step, laser), "&
         &"default=smoothstep", .true.)

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

    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(n_dim, n_cond))
    allocate(ic%seed_r1(n_dim, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    call CFG_get(cfg, "domain_len", dlen)
    ic%seed_r0 = dlen * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = dlen * reshape(tmp_vec, [n_dim, n_cond])

    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "stochastic_density", ic%stochastic_density)
    call CFG_get(cfg, "seed_density", ic%seed_density)
    call CFG_get(cfg, "seed_charge_type", ic%seed_charge_type)
    call CFG_get(cfg, "seed_width", ic%seed_width)
    call CFG_get(cfg, "seed_falloff", ic%seed_falloff)

    init_conds = ic

  end subroutine init_cond_initialize

  !> Add a stochastic background density to the electrons and ions. Note: this
  !> routine temporarily uses variable i_rhs
  subroutine init_cond_stochastic_density(tree)
    use m_a$D_ghostcell, only: a$D_bc_neumann_zero
    use m_a$D_prolong
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
          call a$D_gc_box(tree%boxes, id, i_rhs, &
               a$D_gc_interp, a$D_bc_neumann_zero)
       end do

       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a$D_prolong_linear_from(tree%boxes, id, i_rhs, add=.true.)
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
    call a$D_restrict_tree(tree, i_electron)
    call a$D_restrict_tree(tree, i_pos_ion)
    call a$D_gc_tree(tree, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero)
    call a$D_gc_tree(tree, i_pos_ion, a$D_gc_interp_lim, a$D_bc_neumann_zero)

  end subroutine init_cond_stochastic_density

  subroutine set_stochastic_density(box)
    use omp_lib
    use m_streamer, only: ST_prng
    type(box$D_t), intent(inout) :: box
    integer                      :: proc_id, IJK
    real(dp)                     :: density

    proc_id = 1+omp_get_thread_num()

    do KJI_DO(1,box%n_cell)
       density = ST_prng%rngs(proc_id)%unif_01() * &
            init_conds%stochastic_density
       box%cc(IJK, i_rhs) = density
    end do; CLOSE_DO
  end subroutine set_stochastic_density

  !> Sets the initial condition
  subroutine init_cond_set_box(box)
    use m_geometry
    type(box$D_t), intent(inout) :: box
    integer                     :: IJK, n, nc
    real(dp)                    :: rr($D)
    real(dp)                    :: density

    nc = box%n_cell
    box%cc(DTIMES(:), i_electron) = init_conds%background_density
    box%cc(DTIMES(:), i_pos_ion)  = init_conds%background_density
    box%cc(DTIMES(:), i_phi)      = 0 ! Inital potential set to zero

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

  end subroutine init_cond_set_box

end module m_init_cond_$Dd
