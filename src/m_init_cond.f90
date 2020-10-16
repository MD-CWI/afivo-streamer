#include "../afivo/src/cpp_macros.h"
!> Module to help setting up initial conditions
module m_init_cond
  use m_streamer
  use m_af_all
  use m_chemistry
  use m_types

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
     real(dp), allocatable           :: seed_density2(:)
     integer, allocatable            :: seed_charge_type(:)
     real(dp), allocatable           :: seed_width(:)
     character(string_len), allocatable :: seed_falloff(:)
     integer, allocatable            :: seed1_species(:) !< Custom seed 1 species
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected, public :: init_conds

  public :: init_cond_initialize
  public :: init_cond_set_box
  public :: init_cond_stochastic_density

contains

  ! Set the initial conditions from the configuration
  subroutine init_cond_initialize(tree, cfg)
    use m_config
    type(af_t), intent(in)     :: tree
    type(CFG_t), intent(inout) :: cfg !< Settings

    integer                    :: n, n_cond, varsize, empty_int(0)
    real(dp)                   :: empty_real(0)
    real(dp), allocatable      :: tmp_vec(:)
    character(len=name_len)    :: empty_string(0)
    character(len=name_len), allocatable :: seed_species(:)
    type(initcnd_t)            :: ic

    call CFG_add(cfg, "background_density", 0.0_dp, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "stochastic_density", 0.0_dp, &
         "Stochastic background density (1/m3)")
    call CFG_add(cfg, "seed_density", empty_real, &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed_rel_r0", empty_real, &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", empty_real, &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_charge_type", empty_int, &
         "Type of seed: neutral (0), ions (1) or electrons (-1)", .true.)
    call CFG_add(cfg, "seed_width", empty_real, &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed_falloff", empty_string, &
         "Fall-off type for seed (sigmoid, gaussian, smoothstep, step, laser)", .true.)
    call CFG_add(cfg, "seed1_species", empty_string, &
         "Names of custom species for the first seed", .true.)

    call CFG_get_size(cfg, "seed_density", n_cond)
    ic%n_cond = n_cond

    call CFG_get_size(cfg, "seed_rel_r0", varsize)
    if (varsize /= NDIM * n_cond) &
         stop "seed_rel_r0 variable has incompatible size"

    call CFG_get_size(cfg, "seed_rel_r1", varsize)
    if (varsize /= NDIM * n_cond) &
         stop "seed_rel_r1 variable has incompatible size"

    call CFG_get_size(cfg, "seed_charge_type", varsize)
    if (varsize /= n_cond) &
         stop "seed_charge_type variable has incompatible size"

    call CFG_get_size(cfg, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_width variable has incompatible size"

    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_density2(n_cond))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(NDIM, n_cond))
    allocate(ic%seed_r1(NDIM, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))

    allocate(tmp_vec(NDIM * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    ic%seed_r0 = reshape(tmp_vec, [NDIM, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = reshape(tmp_vec, [NDIM, n_cond])

    do n = 1, n_cond
       ic%seed_r0(:, n) = ic%seed_r0(:, n) * ST_domain_len + ST_domain_origin
       ic%seed_r1(:, n) = ic%seed_r1(:, n) * ST_domain_len + ST_domain_origin
    end do

    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "stochastic_density", ic%stochastic_density)
    call CFG_get(cfg, "seed_density", ic%seed_density)
    call CFG_get(cfg, "seed_charge_type", ic%seed_charge_type)
    call CFG_get(cfg, "seed_width", ic%seed_width)
    call CFG_get(cfg, "seed_falloff", ic%seed_falloff)

    ! Keep density at endpoint the same by default
    call CFG_add(cfg, "seed_density2", ic%seed_density, &
         "Initial density of the seed at other endpoint (1/m3)")
    call CFG_get(cfg, "seed_density2", ic%seed_density2)

    call CFG_get_size(cfg, "seed1_species", varsize)

    if (varsize > 0) then
       allocate(seed_species(varsize))
       allocate(ic%seed1_species(varsize))
       call CFG_get(cfg, "seed1_species", seed_species)
       do n = 1, varsize
          ic%seed1_species(n) = af_find_cc_variable(tree, seed_species(n))
       end do
    end if

    init_conds = ic

  end subroutine init_cond_initialize

  !> Add a stochastic background density to the electrons and ions. Note: this
  !> routine temporarily uses variable i_rhs
  subroutine init_cond_stochastic_density(tree)
    use m_af_ghostcell, only: af_bc_neumann_zero
    use m_af_prolong
    type(af_t), intent(inout) :: tree
    integer                    :: my_lvl, lvl, i, id

    if (init_conds%stochastic_density <= 0.0_dp) return

    ! Determine at which level to create the background density. This is the
    ! highest level that is fully refined
    do my_lvl = 1, tree%highest_lvl
       if (size(tree%lvls(my_lvl)%leaves) > 0) exit
    end do

    ! Use i_rhs to store the stochastic density at this level
    call af_tree_clear_cc(tree, i_rhs)
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
          call af_gc_box(tree, id, [i_rhs])
       end do

       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call af_prolong_linear_from(tree%boxes, id, i_rhs, add=.true.)
       end do

    end do

    ! Finally, add the stochastic density to the electrons and ions
    do lvl = my_lvl, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call af_box_add_cc(tree%boxes(id), i_rhs, i_electron)
          call af_box_add_cc(tree%boxes(id), i_rhs, i_1pos_ion)
       end do
    end do

    ! Restrict and fill ghost cells
    call af_restrict_tree(tree, [i_electron, i_1pos_ion])
    call af_gc_tree(tree, [i_electron, i_1pos_ion])

  end subroutine init_cond_stochastic_density

  subroutine set_stochastic_density(box)
    use omp_lib
    use m_streamer, only: ST_prng
    type(box_t), intent(inout) :: box
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
    use m_gas
    use m_user_methods
    type(box_t), intent(inout) :: box
    integer                    :: IJK, n, nc
    real(dp)                   :: rr(NDIM)
    real(dp)                   :: density

    nc = box%n_cell
    box%cc(DTIMES(:), i_electron) = init_conds%background_density
    box%cc(DTIMES(:), i_1pos_ion) = init_conds%background_density

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])

       if (gas_dynamics) then
          if (associated(user_gas_density)) then
             box%cc(IJK, i_gas_dens) = user_gas_density(box, IJK)
          else
             ! Start with a constant gas number density
             box%cc(IJK, i_gas_dens) = gas_number_density
          end if

          ! Initialize Euler variables: density, momentum, energy
          box%cc(IJK, gas_vars(i_rho)) = box%cc(IJK, i_gas_dens) * &
               gas_molecular_weight
          box%cc(IJK, gas_vars(i_mom)) = 0.0_dp
          box%cc(IJK, gas_vars(i_e)) = &
               gas_pressure * 1e5_dp / (gas_euler_gamma - 1) + &
               0.5_dp * sum(box%cc(IJK, gas_vars(i_mom))**2) / &
               box%cc(IJK, gas_vars(i_rho))
       else if (associated(user_gas_density)) then
          box%cc(IJK, i_gas_dens) = user_gas_density(box, IJK)
       end if

       do n = 1, init_conds%n_cond
          density = GM_density_line(rr, init_conds%seed_r0(:, n), &
               init_conds%seed_r1(:, n), &
               init_conds%seed_density(n), init_conds%seed_density2(n), NDIM, &
               init_conds%seed_width(n), &
               init_conds%seed_falloff(n))

          if (n == 1 .and. allocated(init_conds%seed1_species)) then
             box%cc(IJK, init_conds%seed1_species) = &
                  box%cc(IJK, init_conds%seed1_species) + density
          else
             ! Add electrons and/or ions depending on the seed charge type
             select case (init_conds%seed_charge_type(n))
             case (-1)
                box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + density
             case (0)
                box%cc(IJK, i_1pos_ion) = box%cc(IJK, i_1pos_ion) + density
                box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + density
             case (1)
                box%cc(IJK, i_1pos_ion) = box%cc(IJK, i_1pos_ion) + density
             case default
                error stop "Invalid seed_charge_type"
             end select
          end if
       end do
    end do; CLOSE_DO

  end subroutine init_cond_set_box

end module m_init_cond
